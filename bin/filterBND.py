#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse
from anacore.gtf import loadModel
from anacore.vcf import HeaderFilterAttr
from anacore.region import Region, RegionList
from anacore.fusion import getStrand
from anacore.fusionVcf import AnnotBreakendVCFIO, getBNDInterval


########################################################################
#
# FUNCTIONS
#
########################################################################
class AnnotGetter:
    """Class to get genes regions from annotation file."""

    def __init__(self, filepath):
        """
        Build and return an instance of AnnotGetter.

        :param filepath: Path to the annotations file (format: GTF).
        :type filepath: str
        :return: The new instance.
        :rtype: AnnotGetter
        """
        self.filepath = filepath
        self.model = {}

    def getChr(self, chr):
        """
        Return the genes regions on the specified chr.

        :param chr: The chromosome name.
        :type chr: str
        :return: Genes regions on the specified chr.
        :rtype: anacore.region.RegionList
        """
        if chr not in self.model:
            genes = loadModel(self.filepath, "genes", chr)
            self.model[chr] = genes
        return self.model[chr]


def isIG(first, second, annotation_field):
    """
    Return True if one of the two breakends falls in a gene of immunoglobulin.

    :param first: The breakend of the first shard in fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: The breakend of the second shard in fusion.
    :type second: anacore.vcf.VCFRecord
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :return: True if one of the two breakend falls in a gene of immunoglobulin.
    :rtype: boolean
    """
    is_IG = False
    immunoglobulin = {"IGH", "IGK", "IGL", "IGS"}
    record_gene = {annot["SYMBOL"][0:3] for annot in first.info[annotation_field]}
    mate_gene = {annot["SYMBOL"][0:3] for annot in second.info[annotation_field]}
    if len(immunoglobulin & (record_gene | mate_gene)) > 0:
        is_IG = True
    return is_IG


def isHLA(first, second, annotation_field):
    """
    Return True if one of the two breakends falls in a gene of human leukocyte antigen.

    :param first: The breakend of the first shard in fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: The breakend of the second shard in fusion.
    :type second: anacore.vcf.VCFRecord
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :return: True if one of the two breakend falls in a gene of human leukocyte antigen.
    :rtype: boolean
    """
    is_HLA = False
    record_gene = {annot["SYMBOL"][0:4] for annot in first.info[annotation_field]}
    mate_gene = {annot["SYMBOL"][0:4] for annot in second.info[annotation_field]}
    if "HLA-" in record_gene | mate_gene:
        is_HLA = True
    return is_HLA


def isInner(first, second, annotation_field):
    """
    Return True if the two breakends falls in the same gene (in at least one overlapped gene).

    :param first: The breakend of the first shard in fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: The breakend of the second shard in fusion.
    :type second: anacore.vcf.VCFRecord
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :return: True if the two breakends falls in the same gene (in at least one overlapped gene).
    :rtype: boolean
    """
    is_inner_gene = False
    record_gene = {annot["SYMBOL"] for annot in first.info[annotation_field]}
    mate_gene = {annot["SYMBOL"] for annot in second.info[annotation_field]}
    full_overlapping = record_gene & mate_gene
    if len(full_overlapping) > 0:
        first_strand = getStrand(first, True)
        second_strand = getStrand(second, False)
        if first_strand != second_strand:
            is_inner_gene = True
        else:
            strand = "1" if first_strand == "+" else "-1"
            genes_strands = {annot["STRAND"] for annot in first.info[annotation_field] + second.info[annotation_field] if annot["SYMBOL"] in full_overlapping}
            if strand in genes_strands:
                is_inner_gene = True
    return is_inner_gene


def isReadthrough(up, down, annotation_field, genes, rt_max_dist):
    """
    Return True if the two breakends can be a readthrough.

    :param up: The breakend of the first shard in fusion.
    :type up: anacore.vcf.VCFRecord
    :param down: The breakend of the second shard in fusion.
    :type down: anacore.vcf.VCFRecord
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :param genes: The genes regions by chr.
    :type genes: AnnotGetter
    :param rt_max_dist: Maximum distance to evaluate if the fusion is a readthrough.
    :type rt_max_dist: int
    :return: True if the two breakends can be a readthrough.
    :rtype: boolean
    """
    is_readthrough = False
    if up.chrom == down.chrom:
        up_strand = getStrand(up, True)
        down_strand = getStrand(down, False)
        if (up_strand == "+" and down_strand == "+") or (up_strand == "-" and down_strand == "-"):  # Readthrough are +/+ or -/-
            first = up
            second = down
            if first.pos > second.pos:
                first = down
                second = up
            first_start, first_end = getBNDInterval(first)
            second_start, second_end = getBNDInterval(second)
            interval_start = min(first_start, second_start)
            interval_end = max(first_end, second_end) + 1
            if interval_end - interval_start <= rt_max_dist:
                first_bp_gene = {annot["SYMBOL"] for annot in first.info[annotation_field]}
                second_bp_gene = {annot["SYMBOL"] for annot in second.info[annotation_field]}
                full_overlapping_gene = first_bp_gene & second_bp_gene
                only_first_bp_gene = first_bp_gene - second_bp_gene
                only_second_bp_gene = second_bp_gene - first_bp_gene
                if len(only_first_bp_gene) != 0 and len(only_second_bp_gene) != 0:
                    strand_by_gene = {annot["SYMBOL"]: annot["STRAND"] for annot in first.info[annotation_field] + second.info[annotation_field]}
                    strand = "1" if up_strand == "+" else "-1"
                    only_first_bp_gene = {gene for gene in only_first_bp_gene if strand_by_gene[gene] == strand}
                    only_second_bp_gene = {gene for gene in only_second_bp_gene if strand_by_gene[gene] == strand}
                    possible_on_strand = len(only_first_bp_gene) != 0 and len(only_second_bp_gene) != 0
                    if possible_on_strand:
                        interval_region = Region(interval_start, interval_end, up_strand, first.chrom)
                        overlapped_genes = genes.getChr(first.chrom).getOverlapped(interval_region)
                        overlapped_genes = RegionList([gene for gene in overlapped_genes if gene.name not in full_overlapping_gene and gene.strand == up_strand])
                        overlapped_genes_by_name = {gene.name: gene for gene in overlapped_genes}
                        contradict_readthrough = False
                        for start_gene_name in only_first_bp_gene:
                            start_gene = overlapped_genes_by_name[start_gene_name]
                            for end_gene_name in only_second_bp_gene:
                                end_gene = overlapped_genes_by_name[end_gene_name]
                                for interval_gene in overlapped_genes:
                                    if interval_gene.name != start_gene.name and interval_gene.name != end_gene.name:
                                        if not interval_gene.hasOverlap(start_gene) and not interval_gene.hasOverlap(end_gene):
                                            contradict_readthrough = True
                        is_readthrough = not contradict_readthrough
    return is_readthrough


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter readthrough, inner fusions and fusions occuring on HLA or IG.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_filter = parser.add_argument_group('Filters')  # Filters
    group_filter.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag": a tag is added in FILTER field if the fusion fits a filter. In mode "remove": the fusion is removed from the output if it fits filter. [Default: %(default)s]')
    group_filter.add_argument('-r', '--rt-max-dist', default=50000, type=int, help='Maximum distance to evaluate if the fusion is a readthrough. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the file containing variants annotated with anacore-utils/annotBND.py (format: VCF).')
    group_input.add_argument('-a', '--input-annotations', required=True, help='Path to the genome annotations file used with anacore-utils/annotBND.py (format: GTF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the filtered file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_fusions = 0
    nb_filtered = 0
    genes = AnnotGetter(args.input_annotations)
    with AnnotBreakendVCFIO(args.input_variants, "r", args.annotation_field) as reader:
        with AnnotBreakendVCFIO(args.output_variants, "w") as writer:
            # Header
            writer.copyHeader(reader)
            writer.filter["IG"] = HeaderFilterAttr("IG", "One breakend is located on immunoglobulin.")
            writer.filter["HLA"] = HeaderFilterAttr("HLA", "One breakend is located on HLA.")
            writer.filter["Inner"] = HeaderFilterAttr("Inner", "The two breakends are located in the same gene.")
            writer.filter["Readthrough"] = HeaderFilterAttr("Readthrough", "The fusion is readthrough (it concerns the two following genes in the same strand in an interval <= {}).".format(args.rt_max_dist))
            writer.writeHeader()
            # Records
            for first, second in reader:
                nb_fusions += 1
                new_filters = set()
                # In immunoglobulin
                if isIG(first, second, args.annotation_field):
                    new_filters.add("IG")
                # In HLA
                if isHLA(first, second, args.annotation_field):
                    new_filters.add("HLA")
                # Inner gene
                if isInner(first, second, args.annotation_field):
                    new_filters.add("Inner")
                # # Readthrough
                if isReadthrough(first, second, args.annotation_field, genes, args.rt_max_dist):
                    new_filters.add("Readthrough")
                if len(new_filters) != 0:
                    nb_filtered += 1
                # Write result
                if args.mode == "tag":  # Tag mode
                    for record in (first, second):
                        filters = new_filters
                        if len(record.filter) != 0 and record.filter[0] != "PASS":
                            filters = set(record.filter) | new_filters
                        record.filter = sorted(filters)
                        if len(record.filter) == 0:
                            record.filter = ["PASS"]
                        writer.write(record)
                elif len(new_filters) == 0:  # Filter mode and is tagged
                    for record in (first, second):
                        if len(record.filter) == 0:
                            record.filter = ["PASS"]
                        writer.write(record)

    # Log process
    log.info(
        "{:.2%} of fusions have been {} ({}/{})".format(
            0 if nb_fusions == 0 else nb_filtered / nb_fusions,
            "tagged" if args.mode == "tag" else "removed",
            nb_filtered,
            nb_fusions
        )
    )
    log.info("End of job")
