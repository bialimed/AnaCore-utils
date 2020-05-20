#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import pysam
import logging
import argparse
from statistics import median
from anacore.bed import getSortedAreasByChr


########################################################################
#
# FUNCTIONS
#
########################################################################
def getEndOffset(read, region):
    """
    Return the offset between end of read alignment and region corresponding border.

    :param read: Evaluated read.
    :type read: pysam.AlignedSegment
    :param region: Evaluated region source of the read.
    :type region: anacore.region.Region
    :return: The offset between end of read alignment and region corresponding border.
    :return: int
    """
    offset = None
    if read.is_reverse:
        read_aln_end = read.reference_start + 1  # 0-based
        offset = read_aln_end - region.start
    else:
        read_aln_end = read.reference_end  # 1-based
        offset = read_aln_end - region.end
    return offset


def getStartOffset(read, region):
    """
    Return the offset between start of read alignment and region corresponding border.

    :param read: Evaluated read.
    :type read: pysam.AlignedSegment
    :param region: Evaluated region source of the read.
    :type region: anacore.region.Region
    :return: The offset between end of read alignment and region corresponding border.
    :return: int
    """
    offset = None
    if read.is_reverse:
        read_aln_start = read.reference_end  # 1-based
        offset = read_aln_start - region.end
    else:
        read_aln_start = read.reference_start + 1  # 0-based
        offset = read_aln_start - region.start
    return offset


def getOffsetPenalty(read, region, end_is_usable=False):
    """
    Return the offset penalty for the read coming from the region.

    :param read: Evaluated read.
    :type read: pysam.AlignedSegment
    :param region: Evaluated region source of the read.
    :type region: anacore.region.Region
    :param end_is_usable: True if the offset of the end of read is usable.
    :type end_is_usable: bool
    :return: Offset penalty for the read coming from the region.
    :return: float
    """
    # Start
    start_offset = getStartOffset(read, region)
    start_weight = 1
    if (start_offset > 0 and not read.is_reverse) or (start_offset < 0 and read.is_reverse):
        start_weight = 0.99  # Extend before primer is less plausible than after
    # End
    end_offset = 0
    end_weight = 0.98  # Start carries a little weightier than end of the read in decision
    if end_is_usable:
        end_offset = getEndOffset(read, region)
    return abs(start_offset) * start_weight + abs(end_offset) * end_weight


def endOffsetIsUsable(read, first_region, second_region):
    """
    Return true if the read has sufficient length ad quality to overlap all the two regions. In this case The end offset can be used in penalty calculation for best source selection.

    :param read: Evaluated read.
    :type read: pysam.AlignedSegment
    :param first_region: Region evaluated in best source selection.
    :type first_region: anacore.region.Region
    :param second_region: Region evaluated in best source selection.
    :type second_region: anacore.region.Region
    :return: True if the offset of the end of read is usable in offset penalty calculation.
    :return: bool
    """
    end_is_usable = False
    read_length = read.infer_read_length()
    if read_length is not None:
        longer_region_len = max(first_region.length(), second_region.length())
        if longer_region_len <= read_length:  # Read can overlap the entire region
            if read.is_reverse:
                if read.query_alignment_start == 0:  # The end of the read is aligned
                    end_is_usable = True
                elif read.reference_start + 1 <= min(first_region.start, second_region.start):  # Read alignment ends includes the two regions
                    end_is_usable = True
                else:
                    try:
                        if read.query_qualities:
                            if median(read.query_qualities[max(0, read.query_alignment_start - 5):read.query_alignment_start]) > 22:  # Quality of firsts unaligned is ok (unaligned correspnds to adapter)
                                end_is_usable = True
                    except Exception:
                        pass
            else:  # Read is forward
                if read.query_alignment_end + 1 == read_length:  # The end of the read is aligned
                    end_is_usable = True
                elif read.reference_end >= max(first_region.end, second_region.end):  # Read alignment ends includes the two regions
                    end_is_usable = True
                else:
                    try:
                        if read.query_qualities:
                            if median(read.query_qualities[read.query_alignment_end:read.query_alignment_end + 5]) > 22:  # Quality of firsts unaligned is ok (unaligned correspnds to adapter)
                                end_is_usable = True
                    except Exception:
                        pass
    return end_is_usable


def selectBestSource(read, sources):
    """
    Return the best source for the read from a list of possible sources.

    :param read: Evaluated read.
    :type read: pysam.AlignedSegment
    :param region: Amplion regions that can be source of the read.
    :type region: list
    :return: The best source for the read.
    :return: anacore.region.Region
    """
    selected_region = sources[0]
    read_start_aln = read.reference_end if read.is_reverse else read.reference_start + 1
    for curr_region in sources[1:]:
        skip_penalty = False
        if read.cigartuples is not None and read.cigartuples[0][0] not in {4, 5, 8}:  # Read start nt is aligned (4: soft-clip, 5: hard-clip, 8: difference)
            if read.is_reverse:
                if read_start_aln == curr_region.end and read_start_aln != selected_region.end:  # If read 1 starts strictly on the border of the region
                    skip_penalty = True
                    selected_region = curr_region
            else:
                if read_start_aln == curr_region.start and read_start_aln != selected_region.start:  # If read 1 starts strictly on the border of the region
                    skip_penalty = True
                    selected_region = curr_region
        if not skip_penalty:  # The start of read is ambiguous
            end_is_usable = endOffsetIsUsable(read, selected_region, curr_region)
            selected_penalty = getOffsetPenalty(read, selected_region, end_is_usable)
            curr_penalty = getOffsetPenalty(read, curr_region, end_is_usable)
            if curr_penalty < selected_penalty:
                selected_region = curr_region
    return selected_region


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
    overlapped_regions = list()
    if read.is_reverse:
        read_aln_start = read.reference_end
        for curr_region in regions:
            if read_aln_start < curr_region.start - anchor_offset:
                break
            if read_aln_start <= curr_region.end + anchor_offset:
                if read_aln_start >= curr_region.end - anchor_offset:
                    overlapped_regions.append(curr_region)
    else:
        read_aln_start = read.reference_start + 1
        for curr_region in regions:
            if read_aln_start < curr_region.start - anchor_offset:  # Stop search after read position
                break
            if read_aln_start >= curr_region.start - anchor_offset:
                if read_aln_start <= curr_region.start + anchor_offset:
                    overlapped_regions.append(curr_region)
    selected_region = None
    if len(overlapped_regions) == 1:
        selected_region = overlapped_regions[0]
    elif len(overlapped_regions) > 1:
        selected_region = selectBestSource(read, overlapped_regions)
    return selected_region


def hasOverlapOnZOI(amplicon, first_read, second_read=None, min_cov=10):
    """
    Return True if the read/reads pair overlap the zone of interest.

    :param amplicon: The amplicon region. The start stores the primer start and thickStart store the ZOI start.
    :param amplicon: anacore.region.Region
    :param first_read: The evaluated read in single-end or the first read in paired-end.
    :type first_read: pysam.AlignedSegment
    :param second_read: The evaluated second read in paired-end.
    :type second_read: pysam.AlignedSegment
    :param min_cov: Minimum overlap between read/reads pair and the zone of interest.
    :type min_cov: int
    :return: True if the read/reads pair overlap the zone of interest.
    :rtype: bool
    """
    zoi_cov_len = first_read.get_overlap(amplicon.thickStart, amplicon.thickEnd)
    if second_read is not None:
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


def processSingleReads(aln_reader, panel_regions, RG_id_by_source, args):
    """
    Filter and tag the reads by there amplicon source. This method is designed for single-end reads.

    :param aln_reader: Reader on alignments file.
    :type aln_reader: pysam.AlignmentFile
    :param panel_regions: Amplicons regions by chr.
    :type panel_regions: dict
    :param RG_id_by_source: RG id by amplicon name.
    :type RG_id_by_source: dict
    :param args: The namespace extract from the script arguments.
    :type args: Namespace
    :return: Count of reads by category and by amplicons. Structure: {"eval_order": [...], "count_by_category": {...}, "amplicons_count": {...}}
    :rtype: dict
    """
    # Init log_data
    ct_by_category = {
        "total": 0,
        "unmapped": 0,
        "out_target": 0,
        "invalid_strand": 0,
        "only_primers": 0,
        "valid": 0
    }
    reads_by_RG = {}
    for chr, areas_in_chr in panel_regions.items():
        for curr_area in areas_in_chr:
            reads_by_RG[curr_area.name] = {
                "name": curr_area.name,
                "position": "{}:{}-{}".format(curr_area.chrom, curr_area.start, curr_area.end),
                "count": 0
            }
    # Process
    for curr_read in aln_reader.fetch(until_eof=True):
        if not curr_read.is_secondary and not curr_read.is_supplementary:
            ct_by_category["total"] += 1
            if curr_read.is_unmapped:
                ct_by_category["unmapped"] += 1
            else:
                source_region = None
                if curr_read.reference_name in panel_regions:
                    source_region = getSourceRegion(curr_read, panel_regions[curr_read.reference_name], args.anchor_offset)
                if source_region is None:
                    ct_by_category["out_target"] += 1
                elif args.check_strand and not hasValidStrand(curr_read, source_region):
                    ct_by_category["invalid_strand"] += 1
                else:
                    curr_read.set_tag("RG", RG_id_by_source[source_region.name])
                    if hasOverlapOnZOI(source_region, curr_read, None, args.min_zoi_cov):  # Read overlap ZOI
                        ct_by_category["valid"] += 1
                        reads_by_RG[source_region.name]["count"] += 1
                        FH_out.write(curr_read)
                    else:  # Read is only primer
                        ct_by_category["only_primers"] += 1
    return {
        "eval_order": ["unmapped", "out_target", "invalid_strand", "only_primers", "valid"],
        "count_by_category": ct_by_category,
        "amplicons_count": sorted(
            reads_by_RG.values(),
            key=lambda elt: elt["count"],
            reverse=True
        )
    }


def processPairedReads(aln_reader, panel_regions, RG_id_by_source, args):
    """
    Filter and tag the reads by there amplicon source. This method is designed for paired-end reads.

    :param aln_reader: Reader on alignments file.
    :type aln_reader: pysam.AlignmentFile
    :param panel_regions: Amplicons regions by chr.
    :type panel_regions: dict
    :param RG_id_by_source: RG id by amplicon name.
    :type RG_id_by_source: dict
    :param args: The namespace extract from the script arguments.
    :type args: Namespace
    :return: Count of reads by category and by amplicons. Structure: {"eval_order": [...], "count_by_category": {...}, "amplicons_count": {...}}
    :rtype: dict
    """
    # Init log_data
    ct_by_category = {
        "total": 0,
        "unpaired": 0,
        "pair_unmapped": 0,
        "out_target": 0,
        "invalid_strand": 0,
        "valid_single_read": 0,  # It will be replaced by invalid_pair
        "only_primers": 0,
        "valid": 0
    }
    reads_by_RG = {}
    for chr, areas_in_chr in panel_regions.items():
        for curr_area in areas_in_chr:
            reads_by_RG[curr_area.name] = {
                "name": curr_area.name,
                "position": "{}:{}-{}".format(curr_area.chrom, curr_area.start, curr_area.end),
                "count": 0
            }
    # Process
    valid_reads_by_id = dict()
    for curr_read in aln_reader.fetch(until_eof=True):
        if not curr_read.is_secondary and not curr_read.is_supplementary:
            ct_by_category["total"] += 1
            if not curr_read.is_paired:
                ct_by_category["unpaired"] += 1
            elif curr_read.is_unmapped or curr_read.mate_is_unmapped:
                ct_by_category["pair_unmapped"] += 1
            else:
                source_region = None
                if curr_read.reference_name in panel_regions:
                    source_region = getSourceRegion(curr_read, panel_regions[curr_read.reference_name], args.anchor_offset)
                if source_region is None:
                    ct_by_category["out_target"] += 1
                elif args.check_strand and not hasValidStrand(curr_read, source_region):
                    ct_by_category["invalid_strand"] += 1
                else:
                    ct_by_category["valid_single_read"] += 1
                    curr_read.set_tag("RG", RG_id_by_source[source_region.name])
                    if curr_read.query_name in valid_reads_by_id:  # Pair is valid
                        prev_read = valid_reads_by_id[curr_read.query_name]
                        if prev_read.get_tag("RG") == RG_id_by_source[source_region.name]:
                            if hasOverlapOnZOI(source_region, prev_read, curr_read, args.min_zoi_cov):  # Reads overlap ZOI
                                ct_by_category["valid"] += 2
                                FH_out.write(prev_read)
                                FH_out.write(curr_read)
                                valid_reads_by_id[curr_read.query_name] = None
                                reads_by_RG[source_region.name]["count"] += 2
                            else:  # Reads are only primers
                                ct_by_category["only_primers"] += 2
                                valid_reads_by_id[curr_read.query_name] = None
                    else:
                        valid_reads_by_id[curr_read.query_name] = curr_read
    ct_by_category["invalid_pair"] = ct_by_category["valid_single_read"] - ct_by_category["only_primers"] - ct_by_category["valid"]
    ct_by_category.pop("valid_single_read", None)
    return {
        "eval_order": ["unpaired", "pair_unmapped", "out_target", "invalid_strand", "invalid_pair", "only_primers", "valid"],
        "count_by_category": ct_by_category,
        "amplicons_count": sorted(
            reads_by_RG.values(),
            key=lambda elt: elt["count"],
            reverse=True
        )
    }


def writeTSVSummary(out_path, data):
    """
    Write summary in TSV file. It contains information about the number of reads out off target, reversed and valid.

    :param out_path: Path to the output file.
    :type out_path: str
    :param data: The metrics stored in summary.
    :type data: dict
    """
    with open(args.output_summary, "w") as writer:
        summary = data["count_by_category"]
        # Summary
        writer.write("##Summary\nCategory\tCount\tRatio\n")
        for item in data["eval_order"]:
            writer.write(
                "{}\t{}\t{:5f}\n".format(
                    item.capitalize(),
                    summary[item],
                    summary[item] / summary["total"]
                )
            )
        # Details
        amplicons_ct = data["amplicons_count"]
        writer.write(
            "\n".join([
                "",
                "",
                "##Valid by targeted region",
                "Target\tPosition\tCount\tRatio\n"
            ])
        )
        for curr in amplicons_ct:
            writer.write(
                "{}\t{}\t{}\t{:5f}\n".format(
                    curr["name"],
                    curr["position"],
                    curr["count"],
                    curr["count"] / summary["valid"]
                )
            )


def writeJSONSummary(out_path, data):
    """
    Write summary in JSON file. It contains information about the number of reads out off target, reversed and valid.

    :param out_path: Path to the output file.
    :type out_path: str
    :param data: The metrics stored in summary.
    :type data: dict
    """
    with open(args.output_summary, "w") as FH_summary:
        FH_summary.write(
            json.dumps(data, default=lambda o: o.__dict__, sort_keys=True)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add RG corresponding to the amplicon source. For one reads pair the amplicon is determined from the position of the first match position of the two reads (primers start positions).')
    parser.add_argument('-f', '--summary-format', default='tsv', choices=['json', 'tsv'], help='The summary format. [Default: %(default)s]')
    parser.add_argument('-d', '--check-strand', action='store_true', help='With this option the strand of amplicons is checked.')
    parser.add_argument('-m', '--single-mode', action='store_true', help='Process single-end alignments.')
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

    # Get panel regions
    panel_regions = getSortedAreasByChr(args.input_panel)

    # Filter reads in panel
    log_data = None
    tmp_aln = args.output_aln + "_tmp.bam"
    with pysam.AlignmentFile(args.input_aln, "rb") as FH_in:
        RG_id_by_source = dict()
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
            if args.single_mode:
                log_data = processSingleReads(FH_in, panel_regions, RG_id_by_source, args)
            else:
                log_data = processPairedReads(FH_in, panel_regions, RG_id_by_source, args)

    # Sort output file
    pysam.sort("-o", args.output_aln, tmp_aln)
    pysam.index(args.output_aln)
    os.remove(tmp_aln)

    # Write summary
    if args.output_summary is not None:
        if args.summary_format == "json":
            writeJSONSummary(args.output_summary, log_data)
        else:
            writeTSVSummary(args.output_summary, log_data)
    log.info("End of job")
