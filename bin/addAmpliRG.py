#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.7.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import pysam
import logging
import argparse
from anacore.bed import getAreas


########################################################################
#
# FUNCTIONS
#
########################################################################
def getSelectedAreasByChr(input_panel):
    """
    Return chromosome the list of selected areas from a BED file.

    :param input_panel: Path to the amplicons with their primers (format: BED)
    :type input_panel: str
    :return: By chromosome the list of BED's areas. Each area is represented by an instance of Region.
    :rtype: dict
    """
    selected_areas = getAreas(input_panel)
    selected_areas = sorted(selected_areas, key=lambda x: (x.chrom, x.start, x.end))

    area_by_chr = dict()
    for curr_area in selected_areas:
        chrom = curr_area.chrom
        if chrom not in area_by_chr:
            area_by_chr[chrom] = list()
        area_by_chr[chrom].append(curr_area)

    return(area_by_chr)


def getSourceRegion(read, regions, anchor_offset=0):
    """
    Return the region where the read come from. Returns None if no region corresponds to the read.

    :param read: The evaluated read.
    :type read: pysam.AlignedSegment
    :param regions: Evaluated source regions. Each area must be represented by an instance of Region.
    :type regions: list
    :param anchor_offset: The alignment of the read can start at N nucleotids after the start of the primer. This parameter allows to take account the possible mismatches on the firsts read positions.
    :type anchor_offset: int
    :return: The region where the read come from.
    :return: None/anacore.region.Region
    """
    overlapped_region = None
    ref_start = read.reference_start + 1
    ref_end = read.reference_end
    if read.is_reverse:
        for curr_region in regions:
            if ref_end < curr_region.start - anchor_offset:
                break
            if ref_end <= curr_region.end + anchor_offset:
                if ref_end >= curr_region.end - anchor_offset:
                    overlapped_region = curr_region
                    break
    else:
        for curr_region in regions:
            if ref_start < curr_region.start - anchor_offset:
                break
            if ref_start >= curr_region.start - anchor_offset:
                if ref_start <= curr_region.start + anchor_offset:
                    overlapped_region = curr_region
                    break
    return overlapped_region


def pairHasOverlapOnZOI(amplicon, first_read, second_read, min_cov=10):
    """
    Return the region where the read come from. Returns None if no region corresponds to the read.

    :param amplicon: The amplicon region.
    :type amplicon: anacore.region.Region
    :param first_read: The evaluated read1.
    :type first_read: pysam.AlignedSegment
    :param second_read: The evaluated read2.
    :type second_read: pysam.AlignedSegment
    :param min_cov: The minimum cumulative length of reads pair in zone of interest (R1 on ZOI + R2 on ZOI).
    :type min_cov: int
    :return: The region where the read come from.
    :rtype: bool
    """
    zoi_cov_len = first_read.get_overlap(amplicon.thickStart, amplicon.thickEnd)
    zoi_cov_len += second_read.get_overlap(amplicon.thickStart, amplicon.thickEnd)
    return zoi_cov_len >= min_cov


def hasValidStrand(read, ampl_region):
    """
    Return True if the read is stranded like if it comes from the specified region.

    :param read: The evaluated read.
    :type read: pysam.AlignedSegment
    :param ampl_region: The amplicon region.
    :param ampl_region: anacore.region.Region
    :return: True if the read is stranded like if it comes from the specified region.
    :rtype: bool
    """
    has_valid_strand = False
    if read.is_read1:
        if read.is_reverse:
            if ampl_region.strand == "-":
                has_valid_strand = True
        else:
            if ampl_region.strand == "+":
                has_valid_strand = True
    else:
        if read.is_reverse:
            if ampl_region.strand == "+":
                has_valid_strand = True
        else:
            if ampl_region.strand == "-":
                has_valid_strand = True
    return has_valid_strand


def writeTSVSummary(out_path, data):
    """
    Write summary in TSV file. It contains information about the number of reads out off target, reversed and valid.

    :param out_path: Path to the output file.
    :type out_path: str
    :param data: The metrics stored in summary.
    :type data: dict
    """
    with open(args.output_summary, "w") as FH_summary:
        print(
            "Category\tCount\tRatio",
            "Unpaired\t{}\t{:5f}".format(data["unpaired"], data["unpaired"] / data["total"]),
            "Unmapped\t{}\t{:5f}".format(data["pair_unmapped"], data["pair_unmapped"] / data["total"]),
            "Out_target\t{}\t{:5f}".format(data["out_target"], data["out_target"] / data["total"]),
            "Invalid_strand\t{}\t{:5f}".format(data["invalid_strand"], data["invalid_strand"] / data["total"]),
            "Invalid_pair\t{}\t{:5f}".format(data["invalid_pair"], data["invalid_pair"] / data["total"]),
            "Only_primers\t{}\t{:5f}".format(data["only_primers"], data["only_primers"] / data["total"]),
            "Valid\t{}\t{:5f}".format(data["valid"], data["valid"] / data["total"]),
            sep="\n",
            file=FH_summary
        )


def writeJSONSummary(out_path, data):
    """
    Write summary in JSON file. It contains information about the number of reads out off target, reversed and valid.

    :param out_path: Path to the output file.
    :type out_path: str
    :param data: The metrics stored in summary.
    :type data: dict
    """
    eval_order = ["unpaired", "pair_unmapped", "out_target", "invalid_strand", "invalid_pair", "only_primers", "valid"]
    with open(args.output_summary, "w") as FH_summary:
        FH_summary.write(
            json.dumps({"eval_order": eval_order, "results": data}, default=lambda o: o.__dict__, sort_keys=True)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add RG corresponding to the amplicon source. For a reads pair the amplicon is determined from the position of the first match position of the reads (primers start positions).')
    parser.add_argument('-f', '--summary-format', default='tsv', choices=['json', 'tsv'], help='The summary format. [Default: %(default)s]')
    parser.add_argument('-t', '--check-strand', action='store_true', help='With this option the strand of amplicons is checked.')
    parser.add_argument('-l', '--anchor-offset', type=int, default=4, help='The alignment of the read can start at N nucleotids after the start of the primer. This parameter allows to take account the possible mismatches on the firsts read positions. [Default: %(default)s]')
    parser.add_argument('-z', '--min-zoi-cov', type=int, default=10, help='The minimum cumulative length of reads pair in zone of interest. If the number of nucleotids coming from R1 on ZOI + the number of nucleotids coming from R2 on ZOI is lower than this value the pair is counted in "only_primers". [Default: %(default)s]')
    parser.add_argument('-t', '--RG-tag', default='LB', help='RG tag used to store the area ID. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='The path to the alignments files (format: BAM). This file must be sorted by coordinates.')
    group_input.add_argument('-p', '--input-panel', required=True, help='Path to the list of amplicons with their primers (format: BED). Each area must have an unique ID in the name field and a strand. The thickStart field represent the start of ZOI and the thickEnd the end of ZOI.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-aln', required=True, help='The path to the alignments file (format: BAM).')
    group_output.add_argument('-s', '--output-summary', help='The path to the summary file (format: see --summary-format). It contains information about the number of reads out off target, reversed and valid.')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Count
    total_reads = 0.0
    unpaired_reads = 0
    unmapped_pairs = 0
    out_target_reads = 0
    reverse_reads = 0
    valid_reads = 0
    only_primers_reads_valid_pair = 0
    valid_reads_valid_pair = 0

    # Get panel regions
    panel_regions = getSelectedAreasByChr(args.input_panel)

    # Filter reads in panel
    RG_id_by_source = dict()
    tmp_aln = args.output_aln + "_tmp.bam"
    valid_reads_by_id = dict()
    with pysam.AlignmentFile(args.input_aln, "rb") as FH_in:
        # Replace RG in header
        new_header = FH_in.header.to_dict()
        new_header["RG"] = list()
        RG_idx = 1
        for chrom in sorted(panel_regions):
            for curr_area in panel_regions[chrom]:
                new_header["RG"].append({"ID": str(RG_idx), args.RG_tag: curr_area.name})
                RG_id_by_source[curr_area.name] = str(RG_idx)
                RG_idx += 1
        # Parse reads
        with pysam.AlignmentFile(tmp_aln, "wb", header=new_header) as FH_out:
            for curr_read in FH_in.fetch(until_eof=True):
                if not curr_read.is_secondary and not curr_read.is_supplementary:
                    total_reads += 1
                    if not curr_read.is_paired:
                        unpaired_reads += 1
                    elif curr_read.is_unmapped or curr_read.mate_is_unmapped:
                        unmapped_pairs += 1
                    else:
                        source_region = None
                        if curr_read.reference_name in panel_regions:
                            source_region = getSourceRegion(curr_read, panel_regions[curr_read.reference_name], args.anchor_offset)
                        if source_region is None:
                            out_target_reads += 1
                        elif args.check_strand and not hasValidStrand(curr_read, source_region):
                            reverse_reads += 1
                        else:
                            valid_reads += 1
                            curr_read.set_tag("RG", RG_id_by_source[source_region.name])
                            if curr_read.query_name in valid_reads_by_id:  # Pair is valid
                                prev_read = valid_reads_by_id[curr_read.query_name]
                                if prev_read.get_tag("RG") == RG_id_by_source[source_region.name]:
                                    if pairHasOverlapOnZOI(source_region, prev_read, curr_read, args.min_zoi_cov):  # Reads overlap ZOI
                                        valid_reads_valid_pair += 2
                                        FH_out.write(prev_read)
                                        FH_out.write(curr_read)
                                        valid_reads_by_id[curr_read.query_name] = None
                                    else:  # Reads are only primers
                                        only_primers_reads_valid_pair += 2
                                        valid_reads_by_id[curr_read.query_name] = None
                            else:
                                valid_reads_by_id[curr_read.query_name] = curr_read

    # Sort output file
    pysam.sort("-o", args.output_aln, tmp_aln)
    pysam.index(args.output_aln)
    os.remove(tmp_aln)

    # Write summary
    if args.output_summary is not None:
        data = {
            "total": total_reads,
            "unpaired": unpaired_reads,
            "pair_unmapped": unmapped_pairs,
            "out_target": out_target_reads,
            "invalid_strand": reverse_reads,
            "invalid_pair": valid_reads - only_primers_reads_valid_pair - valid_reads_valid_pair,
            "only_primers": only_primers_reads_valid_pair,
            "valid": valid_reads_valid_pair
        }
        if args.summary_format == "json":
            writeJSONSummary(args.output_summary, data)
        else:
            writeTSVSummary(args.output_summary, data)
    log.info("End of job")
