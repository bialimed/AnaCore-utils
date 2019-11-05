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
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse
from anacore.vcf import VCFIO


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='*****************************.')
    parser.add_argument('-s', '--shared-filters', nargs='*', default=["lowAF", "OOT", "homoP", "popAF", "CSQ", "ANN.COLLOC", "ANN.RNA", "ANN.CSQ", "ANN.popAF"], help='Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-sources', required=True, nargs='+', help='Name of the source in same order of --inputs-variants.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', required=True, nargs='+', help='Path to the variants files coming from different callers (format: VCF). The order determine the which AF and AD are retained: the first caller where it is found in this list.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()
    args.shared_filters = set(args.shared_filters)

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    variants = {}
    for curr_caller, curr_file in zip(args.calling_sources, args.input_variants):
        with VCFIO(curr_file) as FH_in:
            for record in FH_in:
                if record.getName() not in variants:
                    variants[record.getName()] = {}
                variants[record.getName()][curr_caller] = {
                    "merged": "MCO_VAR" in record.info,
                    "sub": set(record.info["MCO_VAR"]) if "MCO_VAR" in record.info else []
                }

    selected_variants = {}
    for variant, info_by_caller in variants.items():
        if len(info_by_caller) > 1:
            contains_merged = False
            for (caller, info) in info_by_caller.items():
                if info["merged"]:
                    contains_merged = True
            if contains_merged:
                print(variant)
                selected_variants[variant] = info_by_caller

    # Read merged caller and update choix du plus pertinent split ou non split
