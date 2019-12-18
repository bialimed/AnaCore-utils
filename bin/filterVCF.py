#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.5.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse
from anacore.filters import filtersFromDict
from anacore.vcf import VCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filters VCF on criteria described in JSON file.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-f', '--input-filters', required=True, help='The path to the filters file (format: JSON).')
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    filters = None
    with open(args.input_filters) as data_file:
        filters = filtersFromDict(json.load(data_file))
    nb_kept = 0
    nb_variants = 0
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                nb_variants += 1
                if filters.eval(record):
                    nb_kept += 1
                    FH_out.write(record)
    # Log process
    log.info(
        "{:.2%} of variants have been removed ({}/{})".format(
            0 if nb_variants == 0 else (nb_variants - nb_kept) / nb_variants,
            (nb_variants - nb_kept),
            nb_variants
        )
    )
    log.info("End of job")
