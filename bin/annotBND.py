#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


import os
import sys
import logging
import argparse
from anacore.gtf import loadModel
from anacore.genomicRegion import Intron, Exon
from anacore.fusion import BreakendVCFIO, getBNDInterval, getStrand
from anacore.annotVcf import AnnotVCFIO, HeaderInfoAttr
from anacore.region import Region, splittedByRef


########################################################################
#
# FUNCTIONS
#
########################################################################
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
                    curr_annot = {
                        "SYMBOL": curr_gene.name,
                        "Gene": curr_gene.annot["id"],
                        "Feature": curr_transcript.annot["id"],
                        "Feature_type": "Transcript",
                        "STRAND": None,
                        "EXON": None,
                        "INTRON": None,
                        "CDS_position": None,
                        "Protein_position": None,
                        "Codon_position": None
                    }
                    if curr_transcript.strand is not None:
                        curr_annot["STRAND"] = curr_transcript.strand
                    subregion, subregion_idx = curr_transcript.getSubFromRefPos(bnd_region.start)
                    if issubclass(subregion.__class__, Intron):  # On intron
                        curr_annot["INTRON"] = "{}/{}".format(
                            subregion_idx,
                            len(curr_transcript.children) - 1
                        )
                    else:  # On exon
                        curr_annot["EXON"] = "{}/{}".format(
                            subregion_idx,
                            len(curr_transcript.children)
                        )
                        if len(curr_transcript.proteins) > 1:
                            log.error(
                                "The management of several proteins for one transcript is not implemented. The transcript {} contains several proteins {}.".format(curr_transcript, curr_transcript.proteins),
                                exec_info=True
                            )
                        if len(curr_transcript.proteins) > 0:
                            curr_annot["CDS_position"] = curr_transcript.proteins[0].getNtPosFromRefPos(bnd_region.start)
                            if curr_annot["CDS_position"] is not None:
                                curr_annot["Protein_position"], curr_annot["Codon_position"] = curr_transcript.proteins[0].getPosOnRegion(bnd_region.start)
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
    is_before_break = is_before_break[0]
    # Set BND_stream
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


# def annotModel(first, second, annotation_field):
#     for first_annot in first.info[annotation_field]:
#         first_annot["model"] = []  # mate_transcript_id:break_type:frame_type
#         for second_annot in second.info[annotation_field]:
#             if first_annot["GENE_SHARD"] == "up" and second_annot["GENE_SHARD"] == "down":
#                 current_model = second_annot["Feature"]


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
    Return retained spot positions for the first and the second breakend when they contain CIPOS. Choice is based on placement on exons boundaries contained in CIPOS interval.

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
    if len(first_exons_sup_by_pos) == 0 and len(second_exons_sup_by_pos) == 0:  # No shard contain an exon at breakend pos
        return (first.pos, second.pos)
    else:
        if len(second_exons_sup_by_pos) == 0:  # Only the 5' shard contains at least one exon at breakend pos: the most supported exon boundary for the first breakend is retained
            selected_pos = getMostSupported(first_exons_sup_by_pos)
            offset = selected_pos - first.pos
            return (selected_pos, second.pos + offset)
        elif len(first_exons_sup_by_pos) == 0:  # Only the 3' shard contains at least one exon at breakend pos: the most supported exon boundary for the second breakend is retained
            selected_pos = getMostSupported(second_exons_sup_by_pos)
            offset = selected_pos - second.pos
            return (first.pos + offset, selected_pos)
        else:  # The two shards contain at least one exon at breakend pos
            first_offsets = {pos - first.pos for pos in first_exons_sup_by_pos}
            second_offsets = {pos - second.pos for pos in second_exons_sup_by_pos}
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
                    offset = pos - first.pos
                    if offset in common:
                        supp_by_offset[offset] += support
                offset = getMostSupported(supp_by_offset)
            else:  # No common offset between two breakends: the most supported exon boundary for the first breakend is retained
                selected_pos = getMostSupported(first_exons_sup_by_pos)
                offset = selected_pos - first.pos
            return (first.pos + offset, second.pos + offset)


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
        # print("fit positions to exons for {}:{} {}:{} ({}/{}) {}from {}".format(
        #     first.chrom,
        #     first.pos,
        #     second.chrom,
        #     second.pos,
        #     getStrand(first, True),
        #     getStrand(second, False),
        #     "IMPRECISE " if "IMPRECISE" in first.info else "",
        #     first.info["SRC"]
        # ))
        # Try to fit positions to exons boundaries
        first_exons_pos = exonsPos(first, genes_by_chr)
        second_exons_pos = exonsPos(second, genes_by_chr)
        # print("", first_exons_pos, second_exons_pos)
        if "IMPRECISE" in first.info:
            first_selected_pos = first.pos if len(first_exons_pos) == 0 else getMostSupported(first_exons_pos)
            second_selected_pos = second.pos if len(second_exons_pos) == 0 else getMostSupported(second_exons_pos)
        else:
            first_selected_pos, second_selected_pos = selectedPos(first, first_exons_pos, second, second_exons_pos)
        # print("", first_selected_pos, second_selected_pos)
        first.info["ANNOT_POS"] = first_selected_pos
        second.info["ANNOT_POS"] = second_selected_pos
    first.info[annotation_field] = getGeneAnnot(first, genes_by_chr)
    annotGeneShard(first, annotation_field)
    second.info[annotation_field] = getGeneAnnot(second, genes_by_chr)
    annotGeneShard(second, annotation_field)


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
    with AnnotVCFIO(args.output_variants, "w") as writer:
        with BreakendVCFIO(args.input_variants) as reader:
            # Header
            writer.copyHeader(reader)
            writer.ANN_titles = ["SYMBOL", "Gene", "Feature", "Feature_type", "STRAND", "EXON", "INTRON", "CDS_position", "Protein_position", "GENE_SHARD"]
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
                writer.write(first)
                writer.write(second)
    log.info("End of job")
