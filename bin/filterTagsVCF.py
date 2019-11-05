#!/usr/bin/env python3
#
# Copyright (C) 2019 IUCT-O
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse
from anacore.filters import filtersFromDict
from anacore.vcf import VCFIO, HeaderFilterAttr


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add filter tags from criteria described in JSON file.')
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

    # Load filters declaration
    filters = None
    with open(args.input_filters) as data_file:
        filters_json = json.load(data_file)
        for curr_filter in filters_json:
            filters.append(filtersFromDict(curr_filter))

    # Process filters
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            for curr_filter in filters:
                FH_out.filter[curr_filter.name] = HeaderFilterAttr(curr_filter.name, curr_filter.description)
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                for curr_filter in filters:
                    if record.filter is None:
                        record.filter = ["PASS"]
                    if not filters.eval(record):
                        if len(record.filter) == 1 and record.filter[0] == "PASS":
                            record.filter = [curr_filter.name]
                        else:
                            record.filter.append(curr_filter.name)
                        record.filter.append(curr_filter.name)
                FH_out.write(record)

    # Log process
    log.info("End of job")
