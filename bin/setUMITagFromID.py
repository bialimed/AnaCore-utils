#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import argparse
import logging
import os
import pysam
import sys


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='For each read get the UMI sequence from the ID and place it in tag.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_RG = parser.add_argument_group('Reads groups')  # Reads group
    group_RG.add_argument('--umi-tag', default="RX", help='This tag will be used to store UMI sequence. [Default: %(default)s]')
    group_RG.add_argument('--umi-separator', default=":", help='Separator used in read ID to separate real ID and UMI sequence (":" with bcl2fastq, "_" with umi_tools, ...). [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='The path to the alignments file (format: BAM).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-aln', required=True, help='The path to the outputted alignments file (format: BAM).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with pysam.AlignmentFile(args.input_aln, "rb", check_sq=False) as reader:
        with pysam.AlignmentFile(args.output_aln, "wb", header=reader.header.to_dict()) as writer:
            for curr_read in reader.fetch(until_eof=True):
                umi_seq = curr_read.query_name.rsplit(args.umi_separator, 1)[1]
                curr_read.set_tag(args.umi_tag, umi_seq.replace("+", "-"))  # fgbio use hyphen as UMI pair separator
                writer.write(curr_read)
    log.info("End of job")
