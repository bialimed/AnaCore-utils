#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.4.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import logging
import argparse
from anacore.gtf import loadModel
from anacore.genomicRegion import Intron, Exon
from anacore.fusion import BreakendVCFIO, getBNDInterval, getStrand
from anacore.vcf import HeaderInfoAttr
from anacore.region import Region, splittedByRef


########################################################################
#
# FUNCTIONS
#
########################################################################
def shardIsBeforeBND(record):
    """
    Return True if the fused shard is before the breakend.

    :param record: The BND record.
    :type record: anacore.vcf.VCFRecord
    :return: True if the fused shard is before the breakend.
    :rtype: boolean
    """
    # Get position relative to the break
    is_before_break = []
    for bnd_idx, alt in enumerate(record.alt):
        if alt.startswith("[") or alt.startswith("]"):
            is_before_break.append(False)
        else:
            is_before_break.append(True)
    if len(set(is_before_break)) > 1:
        record_name = record.id if record.id is not None else record.getName()
        log.error(
            "The breakend {} has several fusion partners with different break's configuration.".format(record_name),
            exec_info=True
        )
    return is_before_break[0]


def getDistBeforeCDSForward(pos, protein):
    """
    Return exonic distance between breakend and CDS for a breakend in 5'UTR of a protein on forward strand.

    :param pos: Position of the breakend.
    :type pos: int
    :param protein: The protein object.
    :type protein: anacore.genomicRegion.Protein
    :return: Exonic distance between breakend and CDS for breakend in 5'UTR.
    :rtype: int
    """
    cds_dist = 0
    if protein.start > pos:
        cds_start = protein.getCDSFromTranscript()[0].start
        exons = protein.transcript.children
        cds_found = False
        curr_exon = exons[0]
        curr_exon_idx = 0
        while not cds_found:
            if curr_exon.end < pos:  # if exon is before breakend
                curr_exon_idx += 1
                curr_exon = exons[curr_exon_idx]
            elif curr_exon.start <= pos:  # if exon contains the breakend
                if curr_exon.end >= cds_start:  # if exon contains the breakend and the CDS
                    cds_dist = cds_start - pos
                    cds_found = True
                else:  # if exon contains the breakend and the CDS is in another exon
                    cds_dist += curr_exon.end - pos + 1
                    curr_exon_idx += 1
                    curr_exon = exons[curr_exon_idx]
            else:  # if exon is after breakend
                if curr_exon.end < cds_start:  # if exon is after breakend and before CDS
                    cds_dist += curr_exon.length()
                    curr_exon_idx += 1
                    curr_exon = exons[curr_exon_idx]
                else:  # if exon is after breakend and contains CDS
                    cds_dist += cds_start - curr_exon.start
                    cds_found = True
    return cds_dist


def getDistBeforeCDSReverse(pos, protein):
    """
    Return exonic distance between breakend and CDS for a breakend in 5'UTR of a protein on reverse strand.

    :param pos: Position of the breakend.
    :type pos: int
    :param protein: The protein object.
    :type protein: anacore.genomicRegion.Protein
    :return: Exonic distance between breakend and CDS for breakend in 5'UTR.
    :rtype: int
    """
    cds_dist = 0
    if protein.end < pos:
        cds_end = protein.getCDSFromTranscript()[0].end
        exons = protein.transcript.children[::-1]
        cds_found = False
        curr_exon = exons[-1]
        curr_exon_idx = len(exons) - 1
        while not cds_found:
            if curr_exon.start > pos:  # if exon is before breakend
                curr_exon_idx -= 1
                curr_exon = exons[curr_exon_idx]
            elif curr_exon.end >= pos:  # if exon contains the breakend
                if curr_exon.start <= cds_end:  # if exon contains the breakend and the CDS
                    cds_dist = pos - cds_end
                    cds_found = True
                else:  # if exon contains the breakend and the CDS is in another exon
                    cds_dist += pos - curr_exon.start + 1
                    curr_exon_idx -= 1
                    curr_exon = exons[curr_exon_idx]
            else:  # if exon is after breakend
                if curr_exon.start > cds_end:  # if exon is after breakend and before CDS
                    cds_dist += curr_exon.length()
                    curr_exon_idx -= 1
                    curr_exon = exons[curr_exon_idx]
                else:  # if exon is after breakend and contains CDS
                    cds_dist += curr_exon.end - cds_end
                    cds_found = True
    return cds_dist


def getGeneAnnot(record, genes_by_chr):
    """
    Return genomic items overlapped by the BND record.

    :param record: The BND record.
    :type record: anacore.vcf.VCFRecord
    :param genes_by_chr: By chromosomes a tree where nodes are genes, transcripts, protein, exons and CDS.
    :type genes_by_chr: dict
    :return: The list of annotations (one annotation by overlapped transcript).
    :rtype: list
    """
    record_strand = getStrand(record)
    shard_before_bnd = shardIsBeforeBND(record)
    bnd_region = Region(record.info["ANNOT_POS"], None, None, record.chrom, record.getName())
    annotations = []
    if record.chrom in genes_by_chr:
        overlapped_genes = genes_by_chr[record.chrom].getOverlapped(bnd_region)
        for curr_gene in overlapped_genes:
            overlapped_transcripts = curr_gene.children.getOverlapped(bnd_region)
            if len(overlapped_transcripts) == 0:
                log.warn("The breakpoint {} is contained by gene {} but by 0 of these transcripts.".format(bnd_region, curr_gene))
            else:
                for curr_transcript in overlapped_transcripts:
                    if len(curr_transcript.proteins) > 1:
                        log.error(
                            "The management of several proteins for one transcript is not implemented. The transcript {} contains several proteins {}.".format(curr_transcript, curr_transcript.proteins),
                            exec_info=True
                        )
                    if curr_transcript.strand is None:
                        log.error(
                            "The transcript {} has no strand.".format(curr_transcript),
                            exec_info=True
                        )
                    curr_annot = {
                        "SYMBOL": curr_gene.name,
                        "Gene": curr_gene.annot["id"],
                        "Feature": curr_transcript.annot["id"],
                        "Feature_type": "Transcript",
                        "STRAND": curr_transcript.strand,
                        "Protein": "" if len(curr_transcript.proteins) == 0 else curr_transcript.proteins[0].annot["id"],
                        "RNA_ELT_TYPE": None,
                        "RNA_ELT_POS": None,
                        "CDS_position": None,
                        "Protein_position": None,
                        "Codon_position": None
                    }
                    # Intron, exon and CDS posiion
                    subregion, subregion_idx = curr_transcript.getSubFromRefPos(bnd_region.start)
                    if issubclass(subregion.__class__, Intron):  # On intron
                        curr_annot["RNA_ELT_TYPE"] = "intron"
                        curr_annot["RNA_ELT_POS"] = "{}/{}".format(
                            subregion_idx,
                            len(curr_transcript.children) - 1
                        )
                        if len(curr_transcript.proteins) > 0 and curr_transcript.strand == record_strand:
                            curr_protein = curr_transcript.proteins[0]
                            # Get CDS on last implicated exon for first shard and first implicated exon on second shard
                            ref_pos = subregion.end + 1
                            if shard_before_bnd:
                                ref_pos = subregion.start - 1
                            curr_annot["CDS_position"] = curr_protein.getNtPosFromRefPos(ref_pos)
                            if curr_annot["CDS_position"] is None or curr_annot["CDS_position"] == 1 or curr_annot["CDS_position"] == curr_protein.length:
                                curr_annot["CDS_position"] = None
                                curr_annot["RNA_ELT_TYPE"] += "&utr"
                                if curr_protein.strand == "+":
                                    curr_annot["RNA_ELT_POS"] += "&" + ("5prim" if curr_protein.start > bnd_region.start else "3prim")
                                    if curr_protein.start > bnd_region.start:
                                        curr_annot["CDS_DIST"] = getDistBeforeCDSForward(bnd_region.start, curr_protein)
                                else:
                                    curr_annot["RNA_ELT_POS"] += "&" + ("5prim" if curr_protein.end < bnd_region.start else "3prim")
                                    if curr_protein.end < bnd_region.start:
                                        curr_annot["CDS_DIST"] = getDistBeforeCDSReverse(bnd_region.start, curr_protein)
                            else:
                                curr_annot["Protein_position"], curr_annot["Codon_position"] = curr_protein.getPosOnRegion(ref_pos)
                    else:  # On exon
                        nb_exon = len(curr_transcript.children)
                        curr_annot["RNA_ELT_TYPE"] = "exon"
                        curr_annot["RNA_ELT_POS"] = "{}/{}".format(subregion_idx, nb_exon)
                        if bnd_region.start == subregion.start:
                            if subregion_idx == 1 and subregion.strand == "+":  # Start of the first exon
                                curr_annot["RNA_ELT_TYPE"] += "&transcriptStart"
                            elif subregion_idx == nb_exon and subregion.strand == "-":   # End of the last exon
                                curr_annot["RNA_ELT_TYPE"] += "&transcriptEnd"
                            else:
                                curr_annot["RNA_ELT_TYPE"] += "&splice" + ("Acceptor" if subregion.strand == "+" else "Donor")
                        elif bnd_region.start == subregion.end:
                            if subregion_idx == 1 and subregion.strand == "-":  # Start of the first exon
                                curr_annot["RNA_ELT_TYPE"] += "&transcriptStart"
                            elif subregion_idx == nb_exon and subregion.strand == "+":   # End of the last exon
                                curr_annot["RNA_ELT_TYPE"] += "&transcriptEnd"
                            else:
                                curr_annot["RNA_ELT_TYPE"] += "&splice" + ("Acceptor" if subregion.strand == "-" else "Donor")
                        if len(curr_transcript.proteins) > 0:
                            curr_protein = curr_transcript.proteins[0]
                            curr_annot["CDS_position"] = curr_protein.getNtPosFromRefPos(bnd_region.start)
                            # UTR
                            if curr_annot["CDS_position"] is None:
                                curr_annot["RNA_ELT_TYPE"] += "&utr"
                                if curr_transcript.proteins[0].strand == "+":
                                    curr_annot["RNA_ELT_POS"] += "&" + ("5prim" if curr_protein.start > bnd_region.start else "3prim")
                                    if curr_protein.start > bnd_region.start:
                                        curr_annot["CDS_DIST"] = getDistBeforeCDSForward(bnd_region.start, curr_protein)
                                else:
                                    curr_annot["RNA_ELT_POS"] += "&" + ("5prim" if curr_protein.end < bnd_region.start else "3prim")
                                    if curr_protein.end < bnd_region.start:
                                        curr_annot["CDS_DIST"] = getDistBeforeCDSReverse(bnd_region.start, curr_protein)
                            # Protein position
                            else:
                                curr_annot["Protein_position"], curr_annot["Codon_position"] = curr_protein.getPosOnRegion(bnd_region.start)
                    # Add to annotations
                    annotations.append(curr_annot)
    return annotations


def annotGeneShard(record, annotation_field):
    """
    Add which shard of genes are in fusion (up or down).

    :param record: The annotated BND record. The BND is previously annotated by genomic regions (see getRegionGeneAnnot()).
    :type record: anacore.vcf.VCFRecord
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    """
    is_before_break = shardIsBeforeBND(record)
    for annot in record.info[annotation_field]:
        annot["GENE_SHARD"] = None
        if annot["STRAND"] is not None:
            annot["GENE_SHARD"] = "down"
            if annot["STRAND"] == "+":
                if is_before_break:
                    annot["GENE_SHARD"] = "up"
            else:
                if not is_before_break:
                    annot["GENE_SHARD"] = "up"


def annotModelRetIntron(first, second, annotation_field):
    """
    Add GENE_SHARD and IN_FRAME in annotations in a context where introns may
    have been retained.
    GENE_SHARD determines if the part of the transcript implicated on the fusion
    RNA is the 5' of the original (up) or the 3' (down).
    IN_FRAME determines if the transcript in the 3' shard of the fusion transcript
    express a part of the original protein: 5' shard imports promoter and start
    of the first transcript in right strand, 3' shard imports the end of the
    second transcript in right strand and the phase of the second transcript is
    kept.
    The following table details IN_FRAME values for all the analysed configurations:
    5' shard     3' shard   Inframe   Note
    5'UTR        5'UTR      1         The first does not start
    CDS          5'UTR      ?/1       1 if first BP is on intron or on end of exon and the second is on intron or start of exon and UTR of second as length compatible to phase.
    3'UTR        5'UTR      ?         If the trancription terminator is cut and the sequence continue to exon and splice to the next or if is readthrough
    5'UTR        CDS        ?         Can be use an other TSS and traduction start.
    CDS          CDS        0/1       Check the phase (end of first shard and start of second shard)
    3'UTR        CDS        ?         If the trancription terminator is cut and the sequence continue to exon and splice to the next or if is readthrough
    *            3'UTR      0         The second is not expressed
    non-coding   5'UTR      1         1 if first and second BND are in intron or in splice site.

    :param first: Breakend of the 5' shard of the fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: Breakend of the 3' shard of the fusion.
    :type second: anacore.vcf.VCFRecord
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    """
    first_strand = getStrand(first, True)
    second_strand = getStrand(second, False)
    annotGeneShard(first, annotation_field)
    annotGeneShard(second, annotation_field)
    for second_annot in second.info[annotation_field]:
        if "IN_FRAME" not in second_annot:
            second_annot["IN_FRAME"] = []
    for first_annot in first.info[annotation_field]:
        if "IN_FRAME" not in first_annot:
            first_annot["IN_FRAME"] = []
        for second_annot in second.info[annotation_field]:
            inframe = "0"
            if first_annot["GENE_SHARD"] == "up" and second_annot["GENE_SHARD"] == "down":
                if first_strand == first_annot["STRAND"] and second_strand == second_annot["STRAND"]:
                    if second_annot["Protein"] != "":  # The second RNA is coding
                        inframe = "."
                        if first_annot["Protein"] != "":  # First and second transcripts are coding
                            if not first_annot["RNA_ELT_TYPE"].endswith("utr") and not second_annot["RNA_ELT_TYPE"].endswith("utr"):  # first: CDS and second: CDS
                                is_intron_jct = "intron" in first_annot["RNA_ELT_TYPE"] and "intron" in second_annot["RNA_ELT_TYPE"]
                                is_exon_jct = None
                                if not is_intron_jct:
                                    is_exon_jct = "exon" in first_annot["RNA_ELT_TYPE"] and "exon" in second_annot["RNA_ELT_TYPE"]
                                if is_intron_jct or is_exon_jct:
                                    # ################## TODO: check stop codon ?
                                    inframe = "0"
                                    if (first_annot["Codon_position"] == 3 and second_annot["Codon_position"] == 1) or (second_annot["Codon_position"] - first_annot["Codon_position"] == 1):
                                        inframe = "1"
                            else:  # At least one breakend falls in UTR
                                if second_annot["RNA_ELT_TYPE"].endswith("utr") and second_annot["RNA_ELT_POS"].endswith("3prim"):  # first: * and second: 3'UTR
                                    inframe = "0"
                                elif first_annot["RNA_ELT_TYPE"].endswith("utr") and second_annot["RNA_ELT_TYPE"].endswith("utr"):  # first: UTR and second: UTR
                                    if first_annot["RNA_ELT_POS"].endswith("5prim") and second_annot["RNA_ELT_POS"].endswith("5prim"):  # first: 5'UTR and second: 5'UTR
                                        inframe = "1"
                                elif not first_annot["RNA_ELT_TYPE"].endswith("utr") and second_annot["RNA_ELT_TYPE"].endswith("utr"):  # first: CDS and second: UTR
                                    if second_annot["RNA_ELT_POS"].endswith("5prim"):  # first: CDS and second: 5'UTR
                                        skip = False
                                        if "exon" in first_annot["RNA_ELT_TYPE"] and "intron" in second_annot["RNA_ELT_TYPE"]:
                                            skip = True
                                            if "spliceDonor" in first_annot["RNA_ELT_TYPE"] or "transcriptEnd" in first_annot["RNA_ELT_TYPE"]:
                                                skip = False
                                        elif "intron" in first_annot["RNA_ELT_TYPE"] and "exon" in second_annot["RNA_ELT_TYPE"]:
                                            skip = True
                                            if "spliceAcceptor" in first_annot["RNA_ELT_TYPE"] or "transcriptStart" in first_annot["RNA_ELT_TYPE"]:
                                                skip = False
                                        if not skip:
                                            inframe = "0"
                                            first_and_utr_codon_pos = first_annot["Codon_position"] + (second_annot["CDS_DIST"] % 3)
                                            # ################## TODO: check stop codon ?
                                            if first_and_utr_codon_pos == 3:
                                                inframe = "1"
                        else:  # First is not coding
                            if second_annot["RNA_ELT_TYPE"].endswith("utr") and second_annot["RNA_ELT_POS"].endswith("3prim"):  # first: * and second: 3'UTR
                                inframe = "0"
                            elif second_annot["RNA_ELT_POS"].endswith("5prim"):  # first: non-coding and second: 5'UTR
                                if first_annot["RNA_ELT_TYPE"].startswith("intron"):
                                    if second_annot["RNA_ELT_TYPE"].startswith("intron"):  # from non-coding first in intron to second 5'UTR in intron
                                        inframe = "1"
                                elif "spliceDonor" in first_annot["RNA_ELT_TYPE"] or "transcriptEnd" in first_annot["RNA_ELT_TYPE"]:
                                    if "spliceAcceptor" in second_annot["RNA_ELT_TYPE"] or "transcriptStart" in second_annot["RNA_ELT_TYPE"]:  # from non-coding first on splice donor to second 5'UTR on splice acceptor
                                        inframe = "1"
            first_annot["IN_FRAME"].append(
                "{}:{}".format(second_annot["Feature"], inframe)
            )
            second_annot["IN_FRAME"].append(
                "{}:{}".format(first_annot["Feature"], inframe)
            )
    for first_annot in first.info[annotation_field]:
        first_annot["IN_FRAME"] = "&".join(first_annot["IN_FRAME"])
    for second_annot in second.info[annotation_field]:
        second_annot["IN_FRAME"] = "&".join(second_annot["IN_FRAME"])


def exonsPos(record, genes_by_chr):
    """
    Return by positions of exons boundaries overlapped by the breakend, the number of alternative transcripts with this exon boundaries.

    :param record: Breakdend record with CIPOS.
    :type record: anacore.vcf.VCFRecord
    :param genes_by_chr: By chromosomes a tree where nodes are genes, transcripts, protein, exons and CDS.
    :type genes_by_chr: dict
    :return: By positions of exons boundaries overlapped by the breakend, the number of alternative transcripts with this exon boundaries.
    :rtype: dict
    """
    record_strand = getStrand(record)
    exons_pos = {}
    start, end = getBNDInterval(record)
    interval_region = Region(start, end, None, record.chrom, record.getName())
    if record.chrom in genes_by_chr:
        overlapped_genes = genes_by_chr[record.chrom].getOverlapped(interval_region)
        for curr_gene in overlapped_genes:
            overlapped_transcripts = curr_gene.children.getOverlapped(interval_region)
            for curr_transcript in overlapped_transcripts:
                for subregion in curr_transcript.children.getOverlapped(interval_region):
                    if record_strand == subregion.strand and issubclass(subregion.__class__, Exon):
                        if interval_region.start <= subregion.start and interval_region.end >= subregion.start:  # Breakend match to exon start
                            if subregion.start not in exons_pos:
                                exons_pos[subregion.start] = 1
                            else:
                                exons_pos[subregion.start] += 1
                        if interval_region.start <= subregion.end and interval_region.end >= subregion.end:
                            if subregion.end not in exons_pos:
                                exons_pos[subregion.end] = 1
                            else:
                                exons_pos[subregion.end] += 1
    return exons_pos


def getMostSupported(exons_sup_by_pos):
    """
    Return the position of the most supported exon boundaries.

    :param exons_sup_by_pos: By positions of exons boundaries overlapped by the breakend, the number of alternative transcripts with this exon boundaries.
    :type exons_sup_by_pos: dict
    :return: The position of the most supported exon boundaries.
    :rtype: int
    """
    selected_pos = None
    max_support = 0
    for curr_pos, curr_support in sorted(exons_sup_by_pos.items()):
        if curr_support > max_support:
            max_support = curr_support
            selected_pos = curr_pos
    return selected_pos


def selectedPos(first, first_exons_sup_by_pos, second, second_exons_sup_by_pos):
    """
    Return retained spot positions for the first and the second breakend when they contain CIPOS (exception for imprecise). Choice is based on placement on exons boundaries contained in CIPOS interval.

    :param first: Breakend of the 5' shard of the fusion.
    :type first: anacore.vcf.VCFRecord
    :param first_exons_sup_by_pos: By positions of exons boundaries overlapped by the first breakend, the number of alternative transcripts with this exon boundaries.
    :type first_exons_sup_by_pos: dict
    :param second: Breakend of the 3' shard of the fusion.
    :type second: anacore.vcf.VCFRecord
    :param second_exons_sup_by_pos: By positions of exons boundaries overlapped by the second breakend, the number of alternative transcripts with this exon boundaries.
    :type second_exons_sup_by_pos: dict
    :return: Retained spot positions for the first and the second breakend when they contain CIPOS.
    :rtype: (int, int)
    """
    first_strand = getStrand(first, True)
    second_strand = getStrand(second, False)
    if len(first_exons_sup_by_pos) == 0 and len(second_exons_sup_by_pos) == 0:  # No shard contain an exon at breakend pos
        return (first.pos, second.pos)
    else:
        if len(second_exons_sup_by_pos) == 0:  # Only the 5' shard contains at least one exon at breakend pos: the most supported exon boundary for the first breakend is retained
            selected_pos = getMostSupported(first_exons_sup_by_pos)
            offset = selected_pos - first.pos
            if first_strand == second_strand:
                return (selected_pos, second.pos + offset)
            else:
                cipos = 0 if "CIPOS" not in second.info else second.info["CIPOS"][1]
                return (selected_pos, second.pos + cipos - offset)
        elif len(first_exons_sup_by_pos) == 0:  # Only the 3' shard contains at least one exon at breakend pos: the most supported exon boundary for the second breakend is retained
            selected_pos = getMostSupported(second_exons_sup_by_pos)
            offset = selected_pos - second.pos
            if first_strand == second_strand:
                return (first.pos + offset, selected_pos)
            else:
                cipos = 0 if "CIPOS" not in first.info else first.info["CIPOS"][1]
                return (first.pos + cipos - offset, selected_pos)
        else:  # The two shards contain at least one exon at breakend pos
            first_offsets = {pos - first.pos for pos in first_exons_sup_by_pos}
            second_offsets = {pos - second.pos for pos in second_exons_sup_by_pos}
            second_cipos = 0 if "CIPOS" not in second.info else second.info["CIPOS"][1]
            if first_strand != second_strand:
                second_offsets = {abs(second_cipos - offset) for offset in second_offsets}
            common = first_offsets & second_offsets
            if len(common) == 1:  # Only one common offset between two breakends
                offset = min(common)
            elif len(common) > 1:  # Several common offsets between two breakends: the offset with the most sum of support (exons first + exons second) is retained
                supp_by_offset = {}
                for pos, support in first_exons_sup_by_pos.items():
                    offset = pos - first.pos
                    if offset in common:
                        supp_by_offset[offset] = support
                for pos, support in second_exons_sup_by_pos.items():
                    offset = pos - second.pos
                    if first_strand != second_strand:
                        offset = abs(second_cipos - offset)
                    if offset in common:
                        supp_by_offset[offset] += support
                offset = getMostSupported(supp_by_offset)
            else:  # No common offset between two breakends: the most supported exon boundary for the first breakend is retained
                selected_pos = getMostSupported(first_exons_sup_by_pos)
                offset = selected_pos - first.pos
            if first_strand == second_strand:
                return (first.pos + offset, second.pos + offset)
            else:
                return (first.pos + offset, second.pos + second_cipos - offset)


def annot(first, second, genes_by_chr, annotation_field):
    """
    Annot breakends by overlapping transcripts. In breakends with CIPOS, the spot position of the annotation is previously determined by search of exons boundaries in interval (This position is stor in ANNOT_POS).

    :param first: Breakend of the 5' shard of the fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: Breakend of the 3' shard of the fusion.
    :type second: anacore.vcf.VCFRecord
    :param genes_by_chr: By chromosomes a tree where nodes are genes, transcripts, protein, exons and CDS.
    :type genes_by_chr: dict
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    """
    first_start, first_end = getBNDInterval(first)
    if first_start == first_end:
        first.info["ANNOT_POS"] = first.pos
        second.info["ANNOT_POS"] = second.pos
    else:
        # Try to fit positions to exons boundaries
        first_exons_pos = exonsPos(first, genes_by_chr)
        second_exons_pos = exonsPos(second, genes_by_chr)
        if "IMPRECISE" in first.info:
            first_selected_pos = first.pos if len(first_exons_pos) == 0 else getMostSupported(first_exons_pos)
            second_selected_pos = second.pos if len(second_exons_pos) == 0 else getMostSupported(second_exons_pos)
        else:
            first_selected_pos, second_selected_pos = selectedPos(first, first_exons_pos, second, second_exons_pos)
        first.info["ANNOT_POS"] = first_selected_pos
        second.info["ANNOT_POS"] = second_selected_pos
    first.info[annotation_field] = getGeneAnnot(first, genes_by_chr)
    second.info[annotation_field] = getGeneAnnot(second, genes_by_chr)
    annotModelRetIntron(first, second, annotation_field)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Annotate BND in a VCF with content of a GTF.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-annotations', required=True, help='Path to the file containing the annotations of genes and transcript for the reference used in variant calling. (format: GTF).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the file containing variants. (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the annotated file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load annotations
    log.info("Load model from {}.".format(args.input_annotations))
    genes = loadModel(args.input_annotations, "genes")
    genes_by_chr = splittedByRef(genes)

    # Annot variants
    log.info("Annot variants in {}.".format(args.input_variants))
    with BreakendVCFIO(args.output_variants, "w", args.annotation_field) as writer:
        with BreakendVCFIO(args.input_variants) as reader:
            # Header
            writer.copyHeader(reader)
            writer.ANN_titles = ["SYMBOL", "Gene", "Feature", "Feature_type", "Protein", "STRAND", "RNA_ELT_TYPE", "RNA_ELT_POS", "CDS_position", "Protein_position", "GENE_SHARD", "IN_FRAME"]
            writer.info[args.annotation_field] = HeaderInfoAttr(
                id=args.annotation_field,
                type="String",
                number=".",
                description="Consequence annotations. Format: " + "|".join(writer.ANN_titles)
            )
            writer.info["ANNOT_POS"] = HeaderInfoAttr(
                id="ANNOT_POS",
                type="Integer",
                number="1",
                description="Breakend position used in annotation. It take into account CIPOS to give priority to a breakend on exon boundaries."
            )
            writer.writeHeader()
            # Records
            for first, second in reader:
                annot(first, second, genes_by_chr, args.annotation_field)
                writer.write(first, second)
    log.info("End of job")
