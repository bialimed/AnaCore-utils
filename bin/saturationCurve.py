#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sv import SVIO
import argparse
import logging
import os
import random
import sys


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Extract saturation data (number of unique elements by sampling step) from count by element.')
    parser.add_argument('-f', '--fields-separator', default="\t", help='Random seed used for sampling. [Default: tabulation]')
    parser.add_argument('-n', '--number-steps', default=50, type=int, help='Number of sampling step to produce saturation curve. [Default: %(default)s]')
    parser.add_argument('-r', '--random-seed', default=0, type=int, help='Random seed used for sampling. [Default: %(default)s]')
    parser.add_argument('-s', '--samples-names', nargs='+', help='Samples names if they are not on first line in "input-count".')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-counts', required=True, help='Path to file containing counts in columns by elements in rows (format: separated values).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-saturation', required=True, help='Path to saturation file (format: separated values).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get samples
    samples = args.samples_names
    if args.samples_names is None:
        with SVIO(args.input_counts, has_title=True) as reader:
            samples = reader.titles[1:]

    # Process
    random.seed(args.random_seed)
    with SVIO(args.output_saturation, "w", separator=args.fields_separator) as writer:
        writer.titles = ["sample", "nb_sampled", "nb_unique"]
        for spl_idx, spl_name in enumerate(samples):
            log.info("Process sample {}".format(spl_name))
            with SVIO(args.input_counts, separator=args.fields_separator, has_title=(args.samples_names is None)) as reader:
                # Get counts
                total_count = 0
                curr_idx = 0
                elements = list()
                for record in reader:
                    if not record[0].startswith("__"):  # Skip metadata in htseq-count output
                        curr_count = int(record[1 + spl_idx])
                        for idx in range(curr_count):
                            elements.append(curr_idx)
                        if curr_count != 0:
                            curr_idx += 1
                            total_count += curr_count
                # Step size
                step_size = int(total_count / args.number_steps)
                # Saturation
                unique_elt = set()
                random.shuffle(elements)
                for curr_idx, curr_elt in enumerate(elements):
                    unique_elt.add(curr_elt)
                    if (curr_idx + 1) % step_size == 0:
                        log.info(
                            "Sample {} step {}/{}".format(
                                spl_name,
                                int((curr_idx + 1) / step_size),
                                args.number_steps
                            )
                        )
                        writer.write([
                            spl_name,
                            curr_idx + 1,
                            len(unique_elt)
                        ])
    log.info("End of job")
