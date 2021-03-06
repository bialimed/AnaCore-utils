#!/usr/bin/env python3

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
from datetime import datetime
from anacore.illumina import getRunFolderInfo


########################################################################
#
# FUNCTIONS
#
########################################################################
class DateTimeEncoder(json.JSONEncoder):
    """JSON encoder for datetime"""
    def default(self, obj):
        encoded = None
        if isinstance(obj, datetime):
            encoded = obj.timestamp()
        else:
            encoded = json.JSONEncodzer.default(self, obj)
        return encoded


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Dump run information coming from several files in run folder.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-r', '--input-run-folder', required=True, help='Path to the run folder.')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-info', help='Path to the outputted info file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with open(args.output_info, "w") as FH_out:
        json.dump(getRunFolderInfo(args.input_run_folder), FH_out, cls=DateTimeEncoder)
    log.info("End of job")
