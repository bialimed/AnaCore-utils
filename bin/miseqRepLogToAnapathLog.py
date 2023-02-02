#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
import argparse
from anacore.illumina.run import CompletedJobInfo


########################################################################
#
# FUNCTIONS
#
########################################################################
def writeLog(log_path, log):
    """
    @summary: Writes log information in txt.
    @param log_path: [str] Path to the log file.
    @param log: [CompletedJobInfo] The workflo information.
    """
    with open(log_path, "w") as FH_log:
        FH_log.write(
            "Workflow={}\n".format(log.workflow_name) + \
            "Version={}\n".format(log.version) + \
            "Parameters={}\n".format(json.dumps(log.parameters)) + \
            "Start_time={}\n".format(log.start_datetime.timestamp()) + \
            "End_time={}\n".format(log.end_datetime.timestamp())
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Converts MiSeq Reporter log file in Anapath log file.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input', required=True, help='The path to MiSeq Reporter log file (format: XML).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output', required=True, help='The path to the outputted log file (format: txt).')
    args = parser.parse_args()

    # Process
    log = CompletedJobInfo(args.input)
    writeLog(args.output, log)
