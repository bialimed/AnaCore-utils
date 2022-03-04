#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.region import consolidated, Region
from anacore.sv import HashedSVIO
import argparse
import json
import logging
from numpy import quantile
import os
import pysam
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def argparse_type_min_base_qual(qual):
    """
    Return value casted to int and check minimum constraint.

    :param qual: Command line argument value for min_base_qual.
    :type qual: str
    :return: Value casted to int and check minimum constraint.
    :rtype: int
    """
    qual = int(qual)
    if qual < 2:
        raise argparse.ArgumentTypeError("Minimum min-base-qual is 2.")
    return qual


def depthsTargets(aln_path, targets, depth_mode, min_base_qual, min_depths, log):
    """
    Add depths distribution and number of nt below depth thresholds in each target.

    :param aln_path: Path to the alignment file (format: SAM/BAM).
    :type aln_path: str
    :param targets: Targets (name, locations, size and metadata). In one target "locations" must not contain overlap.
    :type targets: list
    :param depth_mode: How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once).
    :type depth_mode: str
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param min_depths: Depth thresholds used to count the number of nt below.
    :type min_depths: list
    :param log: Logger of the script.
    :type log: logging.Logger
    """
    nb_selected_regions = len(targets)
    idx_in_part = 1
    with pysam.AlignmentFile(aln_path, "rb") as FH_bam:
        for idx_region, curr_target in enumerate(targets):
            if idx_in_part > nb_selected_regions / 10:
                idx_in_part = 0
                log.info("Processed regions {}/{}.".format(idx_region + 1, nb_selected_regions))
            idx_in_part += 1
            curr_target["depths"] = {
                "distribution": None,
                "under_thresholds": {elt: 0 for elt in min_depths}
            }
            # Get target depths
            target_depths = []
            if depth_mode == "read":
                for curr_sub_region in curr_target["locations"]:
                    target_depths.extend(
                        getReadsDepths(aln_path, curr_sub_region, min_base_qual)
                    )
            else:
                for curr_sub_region in curr_target["locations"]:
                    target_depths.extend(
                        getFragmentsDepths(FH_bam, curr_sub_region, min_base_qual)
                    )
            # Store depths
            for threshold in curr_target["depths"]["under_thresholds"]:
                curr_target["depths"]["under_thresholds"][threshold] = sum([1 for elt in target_depths if elt < threshold])
            curr_target["depths"]["distribution"] = {
                "lower_quartile": quantile(target_depths, 0.25, interpolation="midpoint"),
                "max": max(target_depths),
                "median": quantile(target_depths, 0.5),
                "min": min(target_depths),
                "upper_quartile": quantile(target_depths, 0.75, interpolation="midpoint"),
            }


def getFragmentsDepths(FH_bam, region, min_base_qual=13):
    """
    Return list of depths on region (in fragments).

    :param FH_bam: File handle to the alignments file.
    :type FH_bam: pysam.AlignmentFile
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param excluded_flags: Discard any read that has any of the flags specified in the comma-separated list.
    :type excluded_flags: str
    :return: Depths on region.
    :rtype: list
    """
    depths = list()
    curr_checked = region.start - 1
    for pileupcolumn in FH_bam.pileup(
        region.reference.name,
        region.start - 1,
        region.end,
        max_depth=100000000,
        ignore_overlaps=False,
        ignore_orphans=False,
        min_base_quality=10
    ):
        if pileupcolumn.reference_pos + 1 >= region.start and pileupcolumn.reference_pos + 1 <= region.end:
            # Missing positions
            while curr_checked < pileupcolumn.reference_pos:
                depths.append(0)
                curr_checked += 1
            # Current position
            curr_frag = set()
            for pileupread in pileupcolumn.pileups:
                if pileupcolumn.reference_pos + 1 < region.start or pileupcolumn.reference_pos + 1 > region.end:
                    raise Exception("The reference position {}:{} is out of target {}.".format(region.reference.name, pileupcolumn.reference_pos, region))
                if not pileupread.alignment.is_secondary and not pileupread.alignment.is_duplicate and not pileupread.is_refskip:
                    curr_frag.add(pileupread.alignment.query_name)
            depths.append(len(curr_frag))
            curr_checked = pileupcolumn.reference_pos + 1
    # Missing positions
    while curr_checked < region.end:
        depths.append(0)
        curr_checked += 1
    return depths


# Fast fragments depth but bug following pysam calls with option -s
# def getFragmentsDepths(in_aln, region, min_base_qual=13, excluded_flags="UNMAP,SECONDARY,QCFAIL,DUP"):
#     """
#     Return list of depths on region (in fragments).
#
#     :param in_aln: Path to the alignments file.
#     :type in_aln: str
#     :param region: The evaluated region.
#     :type region: anacore.region.Region
#     :param min_base_qual: Minimum base quality to count this read base on depth.
#     :type min_base_qual: int
#     :param excluded_flags: Discard any read that has any of the flags specified in the comma-separated list.
#     :type excluded_flags: str
#     :return: Depths on region.
#     :rtype: list
#     """
#     out_str = pysam.depth(
#         "-a",
#         "-G", excluded_flags,
#         "-J",
#         "-q", str(min_base_qual),
#         "-r", "{}:{}-{}".format(region.reference.name, region.start, region.end),
#         "-s",
#         in_aln
#     )
#     return [int(elt.split()[2]) for elt in out_str.strip().split("\n")]


def getReadsDepths(in_aln, region, min_base_qual=13, excluded_flags="UNMAP,SECONDARY,QCFAIL,DUP"):
    """
    Return list of depths on region (in reads).

    :param in_aln: Path to the alignments file.
    :type in_aln: str
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param excluded_flags: Discard any read that has any of the flags specified in the comma-separated list.
    :type excluded_flags: str
    :return: Depths on region.
    :rtype: list
    """
    out_str = pysam.depth(
        "-a",
        "-G", excluded_flags,
        "-J",
        "-q", str(min_base_qual),
        "-r", "{}:{}-{}".format(region.reference.name, region.start, region.end),
        in_aln
    )
    return [int(elt.split()[2]) for elt in out_str.strip().split("\n")]


def getTargets(in_targets):
    """
    Return the list of targets (name, locations, size and metadata) from TSV file.

    :param in_targets: Path to targets file (format: TSV). Required columns area
    "name" and "locations". The others will be add in metadata. each line is a
    target corresponding to one or more sub-locations. The locations column
    contains the list of sub-locations (e.g.: "chr1:100-200,chr1500-680").
    :type in_targets: str
    :return: List of targets (name, locations, size and metadata).
    :rtype: list
    """
    targets = []
    with HashedSVIO(in_targets) as reader:
        for record in reader:
            sub_regions = consolidated(
                [Region.fromStr(elt.strip()) for elt in record["locations"].split(",")]
            )
            targets.append({
                "name": record["name"],
                "locations": sub_regions,
                "size": sum([elt.length() for elt in sub_regions]),
                "metadata": {key: val for key, val in record.items() if key not in {"name", "locations"} and val != ""}
            })
    return targets


def writeJSON(out_path, annotated_targets, depth_mode, min_base_qual, min_depths):
    """
    Write targets and their count below depth thresholds in a JSON file.

    :param out_path: Path to the output file.
    :type out_path: str
    :param annotated_targets: List of targets and their count below depth thresholds.
    :type annotated_targets: list
    :param depth_mode: How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once).
    :type depth_mode: str
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param min_depths: Depth thresholds used to count the number of nt below.
    :type min_depths: list
    """
    with open(out_path, "w") as writer:
        for target in annotated_targets:
            target["locations"] = [str(elt) for elt in target["locations"]]  # Region to string
        writer.write(
            json.dumps(
                {
                    "parameters": {
                        "depth_mode": depth_mode,
                        "min_base_qual": min_base_qual,
                        "min_depths": min_depths
                    },
                    "results": annotated_targets
                },
                sort_keys=True
            )
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Write depths distribution and number of nt below depth thresholds for each target.')
    parser.add_argument('-d', '--min-depths', type=int, nargs='+', default=[25, 30], help='Depth thresholds used to count the number of nt below. [Default: %(default)s]')
    parser.add_argument('-m', '--depth-mode', choices=["read", "fragment"], default="fragment", help='How count the depth: by reads (each reads is added independently) or by fragment (the R1 and R2 coming from the same pair are counted only once). [Default: %(default)s]')
    parser.add_argument('-q', '--min-base-qual', type=argparse_type_min_base_qual, default=10, help="Minimum quality to take a read base into account in depth calculation. Must be greater than 1. [Default: %(default)s]")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-a', '--input-aln', required=True, help='Path to the alignments file (format: BAM).')
    group_input.add_argument('-t', '--input-targets', required=True, help='Path to targets file (format: TSV). Required columns area "name" and "locations". The others will be add in metadata. each line is a target corresponding to one or more sub-locations. The locations column contains the list of sub-locations (e.g.: "chr1:100-200,chr1500-680").')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output', required=True, help='Targets and their count below depth thresholds (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load selected regions
    log.info("Load targeted regions.")
    targets = getTargets(args.input_targets)

    # Find shallow areas
    log.info("Find validate depths rate by target.")
    depthsTargets(
        args.input_aln, targets, args.depth_mode, args.min_base_qual,
        args.min_depths, log
    )

    # Write output
    log.info("Write output.")
    writeJSON(args.output, targets, args.depth_mode, args.min_base_qual, args.min_depths)
    log.info("End of job")
