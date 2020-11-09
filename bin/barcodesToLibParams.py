#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import argparse
from csv import DictReader, DictWriter
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
    parser = argparse.ArgumentParser(description='Write lane config file in picard tools IlluminaBasecallsToSam library_params format form a barcodes definition file in picard tools ExtractIlluminaBarcodes barcode_file format.')
    parser.add_argument('-l', '--lane', required=True, help='Lane number.')
    parser.add_argument('-p', '--spl-file-pattern', default="{}_L{}.bam", help='Path pattern for output files in output-config. [Default: %(default)s]')
    parser.add_argument('-s', '--skip-undetermined', action='store_true', help='With this parameter, the Undetermined library is not added to the end of file.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-barcodes', required=True, help='Path to the barcodes file (format: TSV). See picard tools ExtractIlluminaBarcodes barcode_file format.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-config', required=True, help='Path to the output file (format: TSV). See picard tools IlluminaBasecallsToSam library_params format.')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with open(args.output_config, "w") as handle_out:
        with open(args.input_barcodes) as handle_in:
            reader = DictReader(handle_in, delimiter="\t")
            writer = DictWriter(handle_out, delimiter="\t", fieldnames=["OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME", "BARCODE_1", "BARCODE_2"])
            writer.writeheader()
            for row in reader:
                row["library_name"] = row["library_name"].replace("_", "-")
                writer.writerow({
                    "OUTPUT": args.spl_file_pattern.format(row["library_name"], args.lane),
                    "SAMPLE_ALIAS": row["library_name"],
                    "LIBRARY_NAME": row["library_name"],
                    "BARCODE_1": row["barcode_sequence_1"],
                    "BARCODE_2": row["barcode_sequence_2"]
                })
            if not args.skip_undetermined:
                writer.writerow({
                    "OUTPUT": args.spl_file_pattern.format("Undetermined", args.lane),
                    "SAMPLE_ALIAS": "Undetermined",
                    "LIBRARY_NAME": "Undetermined",
                    "BARCODE_1": "N",
                    "BARCODE_2": "N"
                })

    log.info("End of job")
