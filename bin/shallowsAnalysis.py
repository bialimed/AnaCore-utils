#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.3.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import pysam
import logging
import argparse
from anacore.region import Region, RegionList, splittedByRef, iterOverlappedByRegion
from anacore.genomicRegion import Intron
from anacore.bed import getAreas
from anacore.gtf import loadModel
from anacore.vcf import VCFIO
from anacore.gff import GFF3IO, GFF3Record


########################################################################
#
# FUNCTIONS
#
########################################################################
def getTargets(in_aln, in_targets=None):
    """
    Return the list of targeted regions.

    :param in_aln: Path to the alignment file (format: SAM/BAM).
    :type in_aln: str
    :param in_targets: Path to the targeted regions (format: BED). They must not contains any overlap.
    :type in_targets: str
    :return: List of targeted regions.
    :rtype: anacore.region.RegionList
    """
    selected_regions = RegionList()
    if in_targets is None:
        with pysam.AlignmentFile(in_aln, "rb") as FH_bam:
            for ref_info in FH_bam.header["SQ"]:
                selected_regions.append(
                    Region(1, ref_info["LN"], "+", ref_info["SN"], ref_info["SN"])
                )
    else:
        selected_regions = getAreas(in_targets)
        # Check lack of overlap
        selected_regions = sorted(selected_regions, key=lambda x: (x.reference.name, x.start, x.end))
        prev_region = selected_regions[0]
        for curr_region in selected_regions[1:]:
            if curr_region.reference.name == prev_region.reference.name:
                if prev_region.end >= curr_region.start:
                    raise Exception("The regions {} and {} contains an overlap.".format(prev_region, curr_region))
            prev_region = curr_region
    return selected_regions


def addToShallow(curr_chr, curr_pos, prev_opened, shallows):
    """
    Add current position in current shallow frame if they are consecutive else create a shallow area with previous frame and open new shallow frame with current pos.

    :param curr_chr: Name of the current region.
    :type curr_chr: str
    :param curr_pos: The current position with low DP (0-based).
    :type curr_pos: int
    :param prev_opened: The previous shallow frame ({"start": x, "end": y}).
    :type prev_opened: dict
    :param shallows: The list of shallows areas
    :type shallows: anacore.region.RegionList
    """
    if prev_opened["start"] is None:
        prev_opened["start"] = curr_pos
        prev_opened["end"] = curr_pos
    else:
        if prev_opened["end"] == curr_pos - 1:
            prev_opened["end"] = curr_pos
        else:
            shallows.append(
                Region(prev_opened["start"] + 1, prev_opened["end"] + 1, "+", curr_chr)
            )
            prev_opened["start"] = curr_pos
            prev_opened["end"] = curr_pos


def shallowFromAlignment(aln_path, selected_regions, depth_mode, min_depth, log):
    """
    Return the list of shallow regions from the alignment file.

    :param aln_path: Path to the alignment file (format: SAM/BAM).
    :type aln_path: str
    :param selected_regions: Targeted regions. They must not contains any overlap between them.
    :type selected_regions: anacore.region.RegionList
    :param depth_mode: How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once).
    :type depth_mode: str
    :param min_depth: All the locations with a depth under this value are reported in shallows areas.
    :type min_depth: int
    :param log: Logger of the script.
    :type log: logging.Logger
    :return: List of shallow regions.
    :rtype: anacore.region.RegionList
    """
    shallow = RegionList()
    nb_selected_regions = len(selected_regions)
    idx_in_part = 1
    with pysam.AlignmentFile(aln_path, "rb") as FH_bam:
        for idx_region, region in enumerate(selected_regions):
            if idx_in_part > nb_selected_regions / 10:
                idx_in_part = 0
                log.info("Processed regions {}/{}.".format(idx_region + 1, nb_selected_regions))
            idx_in_part += 1
            prev_opened = {"start": None, "end": None}
            curr_checked = region.start - 1
            for pileupcolumn in FH_bam.pileup(region.reference.name, region.start - 1, region.end, max_depth=100000000):
                if pileupcolumn.reference_pos + 1 >= region.start and pileupcolumn.reference_pos + 1 <= region.end:
                    # Missing positions
                    while curr_checked < pileupcolumn.reference_pos:
                        addToShallow(region.reference, curr_checked, prev_opened, shallow)
                        curr_checked += 1
                    # Current position
                    curr_reads_depth = 0
                    curr_frag = set()
                    for pileupread in pileupcolumn.pileups:
                        if pileupcolumn.reference_pos + 1 < region.start or pileupcolumn.reference_pos + 1 > region.end:
                            raise Exception("The reference position {}:{} is out of target {}.".format(region.reference.name, pileupcolumn.reference_pos, region))
                        if not pileupread.alignment.is_secondary and not pileupread.alignment.is_duplicate and not pileupread.is_refskip:
                            curr_reads_depth += 1
                            curr_frag.add(pileupread.alignment.query_name)
                    curr_depth = curr_reads_depth
                    if depth_mode == "fragment":
                        curr_depth = len(curr_frag)
                    if min_depth > curr_depth:
                        addToShallow(region.reference, pileupcolumn.reference_pos, prev_opened, shallow)
                    curr_checked = pileupcolumn.reference_pos + 1
            # Missing positions
            while curr_checked < region.end:
                addToShallow(region.reference, curr_checked, prev_opened, shallow)
                curr_checked += 1
            if prev_opened["start"] is not None:
                shallow.append(
                    Region(prev_opened["start"] + 1, prev_opened["end"] + 1, "+", region.reference)
                )
    return shallow


def variantsRegionFromVCF(vcf_path, min_count=1, symbol="GENE", hgvsc="CDS", hgvsp="AA", count="CNT"):
    """
    Return the region object corresponding to the known variants in a VCF.

    :param vcf_path: Path to the variants file (format: VCF).
    :type vcf_path: str
    :param min_count: Minimum number of samples where the variant is known in the databases to use its information.
    :type min_count: int
    :param symbol: Tag used in VCF.info to store the symbol of the gene.
    :type symbol: str
    :param hgvsc: Tag used in VCF.info to store the HGVSc.
    :type hgvsc: str
    :param hgvsp: Tag used in VCF.info to store the HGVSp.
    :type hgvsp: str
    :param count: Tag used in VCF.info to store the number of database's samples with this variant.
    :type count: str
    :return: List of variants regions.
    :rtype: anacore.region.RegionList
    """
    variants_region = None
    with VCFIO(vcf_path) as FH_in:
        variants_region = [
            Region(
                record.pos,
                record.pos + len(record.ref),
                None,
                record.chrom,
                record.id,
                {
                    "id": record.id,
                    "gene": ("" if symbol not in record.info else record.info[symbol]),
                    "HGVSp": ("" if hgvsp not in record.info else record.info[hgvsp]),
                    "HGVSc": ("" if hgvsc not in record.info else record.info[hgvsc]),
                    "count": (None if count not in record.info else int(record.info[count]))
                }
            ) for record in FH_in if (symbol not in record.info or "_ENST" not in record.info[symbol]) and (count not in record.info or int(record.info[count]) >= min_count)
        ]
    return RegionList(variants_region)


def setVariantsByOverlap(queries, variants):
    """
    Annotate each query by the list of variants overlapping them.

    :param queries: Regions to annotate.
    :type queries: anacore.region.Region
    :param variants: The list of variants where overlapped variants will be searched.
    :type variants: anacore.region.RegionList
    """
    variants_by_chr = splittedByRef(variants)
    queries_by_chr = splittedByRef(queries)
    for chrom, curr_query, overlapped_subjects in iterOverlappedByRegion(queries_by_chr, variants_by_chr):
        curr_query.annot["VAR"] = []
        for sbjct in overlapped_subjects:
            curr_query.annot["VAR"].append(sbjct)


def setTranscriptsAnnotByOverlap(queries, transcripts):
    """
    Annotate each query by the information coming from the transcripts overlapping them.

    :param region: Regions to annotate.
    :type region: anacore.region.Region
    :param transcripts: The list of transcripts where overlapped transcripts will be searched.
    :type transcripts: anacore.region.RegionList
    """
    transcripts_by_chr = splittedByRef(transcripts)
    queries_by_chr = splittedByRef(queries)
    for chrom, curr_query, overlapped_subjects in iterOverlappedByRegion(queries_by_chr, transcripts_by_chr):
        curr_query.annot["ANN"] = getTranscriptsAnnot(curr_query, overlapped_subjects)


def getTranscriptsAnnot(region, transcripts):
    """
    Return for each overlapped transcript the location of start and end of the query region on the transcript and the protein.

    :param region: The query region.
    :type region: anacore.region.Region
    :param transcripts: List of transcripted overlapped by the query region.
    :type transcripts: anacore.region.RegionList
    :return: List of annotations (one by transcript).
    :rtype: list
    """
    annotations = []
    for curr_tr in transcripts:
        curr_annot = {
            "SYMBOL": curr_tr.parent.name,
            "Gene": curr_tr.parent.annot["id"],
            "Feature": curr_tr.annot["id"],
            "Feature_type": "Transcript",
            "STRAND": None,
            "start_EXON": None,
            "start_INTRON": None,
            "start_Protein_position": None,
            "end_EXON": None,
            "end_INTRON": None,
            "end_Protein_position": None
        }
        # Overlap on upstream
        overlap_start = {
            "tr_ref_pos": region.start,
            "tr_sub_idx": None,
            "tr_sub_type": None,
            "prot_pos": None
        }
        if region.start < curr_tr.start:  # The region starts before the transcript and overlap the transcript
            overlap_start["tr_ref_pos"] = curr_tr.start
            overlap_start["tr_sub_type"] = "EXON"
            overlap_start["tr_sub_idx"] = "{}/{}".format(
                (1 if curr_tr.strand != "-" else len(curr_tr.children)),
                len(curr_tr.children)
            )
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region): ##############################
                protein = curr_tr.proteins[0]
                overlap_start["prot_pos"] = (1 if curr_tr.strand != "-" else protein.aaLength())
        else:
            subregion, subregion_idx = curr_tr.getSubFromRefPos(region.start)
            if issubclass(subregion.__class__, Intron):  # The region starts in an intron
                # Get first pos of next exon
                downstream_exon_idx = subregion_idx + 1 - 1  # 0-based
                if curr_tr.strand == "-":
                    downstream_exon_idx = subregion_idx - 1  # 0-based
                overlap_start["tr_ref_pos"] = curr_tr.children[downstream_exon_idx].start
                overlap_start["tr_sub_type"] = "INTRON"
                overlap_start["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children) - 1)
            else:
                overlap_start["tr_sub_type"] = "EXON"
                overlap_start["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children))
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region): ##############################
                protein = curr_tr.proteins[0]
                if overlap_start["tr_ref_pos"] < protein.start:  # The region overlap an UTR
                    overlap_start["prot_pos"] = (1 if curr_tr.strand != "-" else protein.aaLength())
                else:
                    overlap_start["prot_pos"] = protein.getPosOnRegion(overlap_start["tr_ref_pos"])[0]
        # Overlap on downstream
        overlap_end = {
            "tr_ref_pos": region.end,
            "tr_sub_idx": None,
            "tr_sub_type": None,
            "prot_pos": None
        }
        if region.end > curr_tr.end:  # The region ends after the transcript and overlap the transcript
            overlap_end["tr_ref_pos"] = curr_tr.end
            overlap_end["tr_sub_type"] = "EXON"
            overlap_end["tr_sub_idx"] = "{}/{}".format(
                (len(curr_tr.children) if curr_tr.strand != "-" else 1),
                len(curr_tr.children)
            )
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region): ##############################
                protein = curr_tr.proteins[0]
                overlap_end["prot_pos"] = (protein.aaLength() if curr_tr.strand != "-" else 1)
        else:
            subregion, subregion_idx = curr_tr.getSubFromRefPos(region.end)
            if issubclass(subregion.__class__, Intron):  # The region ends in an intron
                # Get last pos of previous exon
                upstream_exon_idx = subregion_idx - 1  # 0-based
                if curr_tr.strand == "-":
                    upstream_exon_idx = subregion_idx + 1 - 1  # 0-based
                overlap_end["tr_ref_pos"] = curr_tr.children[upstream_exon_idx].end
                overlap_end["tr_sub_type"] = "INTRON"
                overlap_end["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children) - 1)
            else:
                overlap_end["tr_sub_type"] = "EXON"
                overlap_end["tr_sub_idx"] = "{}/{}".format(subregion_idx, len(curr_tr.children))
            # Annot protein
            if len(curr_tr.proteins) > 0 and curr_tr.proteins[0].hasOverlap(region): ##############################
                protein = curr_tr.proteins[0]
                if overlap_end["tr_ref_pos"] > protein.end:  # The region overlap an UTR
                    overlap_end["prot_pos"] = (protein.aaLength() if curr_tr.strand != "-" else 1)
                else:
                    overlap_end["prot_pos"] = protein.getPosOnRegion(overlap_end["tr_ref_pos"])[0]
        # Store info in annotations
        start = overlap_start
        end = overlap_end
        curr_annot["STRAND"] = "1"
        if curr_tr.strand == "-":
            curr_annot["STRAND"] = "-1"
            start = overlap_end
            end = overlap_start
        curr_annot["start_" + start["tr_sub_type"]] = start["tr_sub_idx"]
        curr_annot["start_Protein_position"] = start["prot_pos"]
        curr_annot["end_" + end["tr_sub_type"]] = end["tr_sub_idx"]
        curr_annot["end_Protein_position"] = end["prot_pos"]
        annotations.append(curr_annot)
    return annotations


def region2dict(region, args):
    """
    Return the conversion of region object to dict.

    :param region: The region to convert.
    :type region: anacore.region.Region
    :param args: Parameters used in analysis.
    :type args: NameSpace
    :return: The dict representation of the region.
    :rtype: dict
    """
    return {
        "start": region.start,
        "end": region.end,
        "strand": region.strand,
        "reference": (None if region.reference is None else region.reference.name),
        "name": region.name,
        "annotations": None if args.input_annotations is None else region.annot["ANN"],
        "known_variants": None if len(args.inputs_variants) == 0 else [curr_var.annot for curr_var in region.annot["VAR"]]
    }


def writeJSON(out_path, shallow, args):
    """
    Write shallow areas there annotations in a JSON file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param shallow: The list of shallow areas.
    :type shallow: anacore.region.RegionList
    :param args: Parameters used in analysis.
    :type args: NameSpace
    """
    nb_shallows = len(shallow)
    with open(out_path, "w") as FH_out:
        FH_out.write("{\n")
        FH_out.write(
            '  "parameters": {{"depth_mode": "{}", "min_depth": {}, "use_annotations": {}, "use_variants": {}, "known_min_count": {}}},\n'.format(
                args.depth_mode,
                args.min_depth,
                str(args.input_annotations is not None).lower(),
                str(len(args.inputs_variants) != 0).lower(),
                args.known_min_count
            )
        )
        FH_out.write('  "results": [\n')
        for idx_shallow, curr_shallow in enumerate(shallow):
            curr_json = json.dumps(region2dict(curr_shallow, args), sort_keys=True)
            suffix = "," if idx_shallow + 1 < nb_shallows else ""
            FH_out.write("  {}{}\n".format(curr_json, suffix))
        FH_out.write('  ]\n')
        FH_out.write('}')


def writeGenesJSON(out_path, shallow, args):
    """
    Write the genes, their parts and known variants affected by shallow areas in a JSON file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param shallow: The list of shallow areas.
    :type shallow: anacore.region.RegionList
    :param args: Parameters used in analysis.
    :type args: NameSpace
    """
    genes_data = {}
    # Annotations
    if args.input_annotations is not None:
        for curr_shallow in sorted(shallow, key=lambda x: (x.reference.name, x.start, x.end)):
            for annot in curr_shallow.annot["ANN"]:
                gene_name = annot["SYMBOL"]
                transcript_id = annot["Feature"].split(".")[0]
                if gene_name not in genes_data:
                    genes_data[gene_name] = {"tr": dict(), "var": list()}
                if transcript_id not in genes_data[gene_name]["tr"]:
                    genes_data[gene_name]["tr"][transcript_id] = list()
                stored_tr = genes_data[gene_name]["tr"][transcript_id]
                start = {"type": "exon", "pos": annot["start_EXON"]}
                if start["pos"] is None:
                    start = {"type": "intron", "pos": annot["start_INTRON"]}
                end = {"type": "exon", "pos": annot["end_EXON"]}
                if end["pos"] is None:
                    end = {"type": "intron", "pos": annot["end_INTRON"]}
                if start["type"] == end["type"] and start["pos"] == end["pos"]:
                    msg = "on {} {}".format(start["type"], start["pos"])
                    if len(stored_tr) == 0 or msg != stored_tr[-1]:  # Prevent repeated location (e.g. 2 shallows on same exon)
                        stored_tr.append(msg)
                else:
                    stored_tr.append(
                        "from {} {} to {} {}".format(
                            start["type"], start["pos"], end["type"], end["pos"]
                        )
                    )
    # Variants
    if len(args.inputs_variants) != 0:
        for curr_shallow in shallow:
            for curr_var in sorted(curr_shallow.annot["VAR"], key=lambda elt: elt.annot["count"], reverse=True):
                curr_var = curr_var.annot
                gene_name = curr_var["gene"]
                if gene_name != "" and curr_var["HGVSc"] != "":
                    if gene_name not in genes_data:
                        genes_data[gene_name] = {"tr": dict(), "var": list()}
                    genes_data[gene_name]["var"].append({
                        "HGVS": curr_var["HGVSc"] if curr_var["HGVSp"] is None or curr_var["HGVSp"] == "p.?" else curr_var["HGVSp"],
                        "id": curr_var["id"],
                        "count": curr_var["count"]
                    })
    # write
    with open(out_path, "w") as FH_out:
        FH_out.write(
            json.dumps({
                "parameters": {
                    "depth_mode": args.depth_mode,
                    "min_depth": args.min_depth,
                    "use_annotations": args.input_annotations is not None,
                    "use_variants": len(args.inputs_variants) != 0,
                    "known_min_count": args.known_min_count
                },
                "results": genes_data
            })
        )


def writeGFF(out_path, shallow, args):
    """
    Write shallow areas there annotations in a GFF3 file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param shallow: The list of shallow areas.
    :type shallow: anacore.region.RegionList
    :param args: Parameters used in analysis.
    :type args: NameSpace
    """
    with GFF3IO(out_path, "w") as FH_out:
        for curr_shallow in sorted(shallow, key=lambda x: (x.reference.name, x.start, x.end)):
            record = GFF3Record(
                curr_shallow.reference.name,
                "shallowsAnalysis",
                "experimental_feature",
                curr_shallow.start,
                curr_shallow.end
            )
            if args.input_annotations is not None:
                for idx, annot in enumerate(curr_shallow.annot["ANN"]):
                    fields = []
                    for k, v in sorted(annot.items()):
                        fields.append("{}:{}".format(k, v))
                    record.annot["ann_{}".format(idx + 1)] = "|".join(fields)
            if len(args.inputs_variants) > 0:
                for idx, var_region in enumerate(curr_shallow.annot["VAR"]):
                    fields = []
                    for k, v in sorted(var_region.annot.items()):
                        fields.append("{}:{}".format(k, v))
                    record.annot["var_{}".format(idx + 1)] = "|".join(fields)
            FH_out.write(record)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Extract shallow areas from the alignment are annotate them with genomic features and known variants.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('-m', '--depth-mode', choices=["read", "fragment"], default="fragment", help='How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once). [Default: %(default)s]')
    parser.add_argument('-d', '--min-depth', type=int, default=30, help='All the locations with a depth under this value are reported in shallows areas. [Default: %(default)s]')
    group_known = parser.add_argument_group('Known variants')
    group_known.add_argument('-n', '--known-count-field', default="CNT", help="Field used in known variants database to store the number of database's samples with this variant. [Default: %(default)s]")
    group_known.add_argument('-i', '--known-hgvsc-field', default="CDS", help='Field used in known variants databases to store the HGVSp. [Default: %(default)s]')
    group_known.add_argument('-p', '--known-hgvsp-field', default="AA", help='Field used in known variants databases to store the HGVSp. [Default: %(default)s]')
    group_known.add_argument('-c', '--known-min-count', type=int, default=3, help='Minimum number of samples where the variant is known in the databases to use its information. [Default: %(default)s]')
    group_known.add_argument('-y', '--known-symbol-field', default="GENE", help='Field used in known variants databases to store the symbol of the gene. This field is required in VCF when you generate the output by gene: --output-genes. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-b', '--input-aln', required=True, help='Path to the alignments file (format: BAM).')
    group_input.add_argument('-t', '--input-targets', help='Path to the targeted regions (format: BED). They must not contains any overlap. [Default: all positions defined in the alignment file header]')
    group_input.add_argument('-a', '--input-annotations', help='Path to the file defining transcripts, genes and proteins locations (format: GTF). This file allow to annotate locations on genes and proteins located on shallows areas. [Default: The shallows areas are not annotated]')
    group_input.add_argument('-s', '--inputs-variants', nargs="+", default=[], help='Path(es) to the file(s) defining known variants (format: VCF). This file allow to annotate variant potentially masked because they are on shallows areas. [Default: The variants on shallows areas are not reported]')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-shallow', default="shallow_areas.gff3", help='Path to the file containing shallow areas and there annotations. (format: GFF3 or JSON if file name ends with ".json"). [Default: %(default)s]')
    group_output.add_argument('-g', '--output-genes', help="Path to the file containing genes's shallow areas (format: JSON). For use this option with the option --inputs-variants the databases must provide the gene's symbol of each variant (see --known-symbol-field).")
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load selected regions
    log.info("Load targeted regions.")
    selected_regions = getTargets(args.input_aln, args.input_targets)

    # Find shallow areas
    log.info("Find shallow areas.")
    shallow = shallowFromAlignment(args.input_aln, selected_regions, args.depth_mode, args.min_depth, log)

    # Annotate shallow areas
    if args.input_annotations is not None:
        log.info("Load annotations from {}.".format(args.input_annotations))
        transcripts = loadModel(args.input_annotations, "transcripts")
        log.info("Annotate shallow areas.")
        setTranscriptsAnnotByOverlap(shallow, transcripts)

    # Retrieved known variants potentialy masked in shallow areas
    for curr_input in args.inputs_variants:
        log.info("Load variants from {}.".format(curr_input))
        variant_regions = variantsRegionFromVCF(
            curr_input,
            args.known_min_count,
            args.known_symbol_field,
            args.known_hgvsc_field,
            args.known_hgvsp_field,
            args.known_count_field
        )
        log.info("List potentialy masked mutations.")
        setVariantsByOverlap(shallow, variant_regions)

    # Write output
    log.info("Write output.")
    if args.output_shallow.endswith(".json"):
        writeJSON(args.output_shallow, shallow, args)
    else:
        writeGFF(args.output_shallow, shallow, args)
    if args.output_genes is not None:
        writeGenesJSON(args.output_genes, shallow, args)
    log.info("End of job")
