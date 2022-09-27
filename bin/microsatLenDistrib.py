#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.bed import getAreas
from anacore.msi.locus import Locus, LocusDataDistrib, LocusRes
from anacore.msi.reportIO import ReportIO
from anacore.msi.sample import MSISample
import argparse
import logging
import os
import pysam
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def getMicrosatLengths(aln_reader, microsat, padding=2, keep_duplicates=False):
    """
    Return length distribution for the microsatellite from reads.

    :param aln_reader: Handle to alignments file.
    :type aln_reader: pysam.AlignmentFile
    :param microsat: Microsatellite region.
    :type microsat: anacore.region.Region
    :param padding: Minimum number of nucleotids aligned on each side around the microsatellite to take read into account. This parameter is used to skip reads containing an incomplete microsatellite.
    :type padding: int
    :param keep_duplicates: If True duplicates reads are taking into account in lengths distribution.
    :type keep_duplicates: boolean
    :return: Length distribution for the microsatellite from reads.
    :rtype: dict
    """
    ms_start_0 = microsat.start - 1
    ms_end_0 = microsat.end - 1
    nb_by_length = {}
    for read in aln_reader.fetch(microsat.reference.name, ms_start_0, microsat.end):  # For each read ovelapping current microsat
        if (keep_duplicates or not read.is_duplicate) and not read.is_secondary and not read.is_supplementary:
            if read.reference_start + 1 <= microsat.start - padding and read.reference_end >= microsat.end + padding:  # Align on microsat and padding
                read_ms_length = getReadRegionLength(read, ms_start_0, ms_end_0)
                if read_ms_length in nb_by_length:
                    nb_by_length[read_ms_length] += 1
                else:
                    nb_by_length[read_ms_length] = 1
    return nb_by_length


def getMicrosatStitchedLengths(aln_reader, microsat, padding=2, keep_duplicates=False):
    """
    Return length distribution for the microsatellite from mate consensus reads.

    :param aln_reader: Handle to alignments file.
    :type aln_reader: pysam.AlignmentFile
    :param microsat: Microsatellite region.
    :type microsat: anacore.region.Region
    :param padding: Minimum number of nucleotids aligned on each side around the microsatellite to take read into account. This parameter is used to skip reads containing an incomplete microsatellite.
    :type padding: int
    :param keep_duplicates: If True duplicates reads are taking into account in lengths distribution.
    :type keep_duplicates: boolean
    :return: Length distribution for the microsatellite from mate consensus reads.
    :rtype: dict
    """
    ms_start_0 = microsat.start - 1
    ms_end_0 = microsat.end - 1
    nb_by_length = dict()
    mate = dict()
    for read in aln_reader.fetch(microsat.reference.name, ms_start_0, microsat.end):  # For each read ovelapping current microsat
        if (keep_duplicates or not read.is_duplicate) and not read.is_secondary and not read.is_supplementary:
            if read.reference_start + 1 <= microsat.start - padding and read.reference_end >= microsat.end + padding:  # Align on microsat and padding
                read_ms_length = getReadRegionLength(read, ms_start_0, ms_end_0)
                if read.query_name not in mate:
                    mate[read.query_name] = read_ms_length
                else:
                    if read_ms_length == mate[read.query_name]:
                        if read_ms_length in nb_by_length:
                            nb_by_length[read_ms_length] += 1
                        else:
                            nb_by_length[read_ms_length] = 1
    return nb_by_length


def getReadRegionLength(read, start, end):
    """
    Return read length on the reference region (take into account deletion and
    insertion).

    :param read: Read.
    :type read: pysam.PileupRead
    :param start: Region start 0-based.
    :type start: int
    :param end: Region end 0-based.
    :type end: int
    :return: Length of the region on read.
    :rtype: int
    """
    read_region_len = 0
    none_before_start = 0
    prev_pos = -1
    for pos in read.get_reference_positions(full_length=True):  # Get length of microsat in read
        if pos is None:
            if prev_pos >= start and prev_pos <= end:
                read_region_len += 1
                none_before_start = 0
            else:
                none_before_start += 1
        else:
            if pos >= start and pos <= end:
                read_region_len += 1 + none_before_start
            none_before_start = 0
            prev_pos = pos
    return read_region_len


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Retrieves the reads length distribution for loci.')
    parser.add_argument('-d', '--keep-duplicates', action='store_true', help='Duplicates reads are taking into account in lengths distribution.')
    parser.add_argument('-m', '--method-name', default="aln", help='Name of the method used to produce lengths data. [Default: %(default)s]')
    parser.add_argument('-p', '--padding', type=int, default=2, help='Minimum number of nucleotids aligned on each side around the microsatellite to take read into account. This parameter is used to skip reads containing an incomplete microsatellite. [Default: %(default)s]')
    parser.add_argument('-n', '--sample-name', help='Name of the sample. [Default: alignments file name without extension]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_stitching = parser.add_argument_group('Stitching')
    stitching_exclusion = group_stitching.add_mutually_exclusive_group()
    stitching_exclusion.add_argument('-e', '--reads-stitched', action='store_true', help='Reads has been stitched before alignment.')
    stitching_exclusion.add_argument('-s', '--stitch-count', action='store_true', help='Reads are not stitched before alignment but count must be processed on mate consensus.')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-a', '--input-alignments', required=True, help='Path to the file containing reads or fragment alignment. (format: BAM).')
    group_input.add_argument('-i', '--input-microsatellites', required=True, help='Path to the file containing microsatellites regions. (format: BED).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-results', required=True, help='Path to the results file. (format: JSON).')
    args = parser.parse_args()
    sample_name = args.sample_name if args.sample_name is not None else os.path.basename(args.input_alignments).split(".")[0]

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load microsatellites regions
    microsat = getAreas(args.input_microsatellites)

    # Load microsatellites size distribution
    locus_by_id = dict()
    with pysam.AlignmentFile(args.input_alignments, "rb") as aln_reader:
        for curr_ms in microsat:
            if args.stitch_count:
                nb_by_length = getMicrosatStitchedLengths(aln_reader, curr_ms, args.padding, args.keep_duplicates)
            else:
                nb_by_length = getMicrosatLengths(aln_reader, curr_ms, args.padding, args.keep_duplicates)
            # Create locus entry
            locus_id = "{}:{}-{}".format(curr_ms.reference.name, curr_ms.start - 1, curr_ms.end)
            locus_by_id[locus_id] = Locus(
                locus_id,
                curr_ms.name,
                {
                    args.method_name: LocusRes(
                        status=None,
                        data={
                            "lengths": LocusDataDistrib(
                                nb_by_length,
                                "fragments" if args.reads_stitched or args.stitch_count else "reads"
                            )
                        }
                    )
                }
            )

    # Write results
    ReportIO.write(
        [MSISample(sample_name, locus_by_id)],
        args.output_results
    )
    log.info("End of job")
