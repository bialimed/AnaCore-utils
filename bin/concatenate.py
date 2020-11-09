#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.abstractFile import isGzip
import argparse
import gzip
import logging
import os
import sys


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Concatenate text files.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--inputs', required=True, nargs='+', help='Pathes to files.')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output', required=True, help='Path to the concatenated file.')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    out_open_fct = open
    out_mode = "w"
    if args.output.endswith('.gz'):
        out_open_fct = gzip.open
        out_mode = "wt"
    with out_open_fct(args.output, out_mode) as writer:
        last_line = "\n"
        for curr_in_file in args.inputs:
            if not last_line.endswith("\n"):
                writer.write("\n")  # Start new line for a new file
            in_open_fct = open
            in_mode = "r"
            if isGzip(curr_in_file):
                in_open_fct = gzip.open
                in_mode = "rt"
            with in_open_fct(curr_in_file, in_mode) as reader:
                for line in reader:
                    writer.write(line)
                    last_line = line
    log.info("End of job")
