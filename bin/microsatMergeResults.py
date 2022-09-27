#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.msi.reportIO import ReportIO
import argparse
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def process(args, log):
    """
    Merge multiple MSI ReportIO from the same samples and same loci.

    :param args: The namespace extracted from the script arguments.
    :type args: Namespace
    :param log: The logger of the script.
    :type log: loggin.Logger
    """
    final_report = ReportIO.parse(args.inputs_reports[0])
    final_loci = final_report[0].loci.keys()
    for spl in final_report:
        curr_loci = spl.loci.keys()
        if final_loci != curr_loci:
            raise Exception("Samples from report {} does not contain the same loci.".format(args.inputs_reports[0]))
    for curr_report_path in args.inputs_reports[1:]:
        curr_report = ReportIO.parse(curr_report_path)
        for final_spl, curr_spl in zip(final_report, curr_report):
            if final_spl.name != curr_spl.name:
                raise Exception("Reports from {} does not contain the same samples.".format(args.inputs_reports))
            curr_loci = set(curr_spl.loci.keys())
            if final_loci != curr_loci:
                raise Exception("Reports from {} does not contain the same loci.".format(args.inputs_reports))
            # Sample results
            for method, res in curr_spl.results.items():
                if method in final_spl.results:
                    log.warnings("Results from the method {} exist in multiple reports and will be overwritten by the last occurence. Report files: {}.".format(method, args.inputs_reports))
                final_spl.results[method] = res
            # Loci results
            for locus_id, locus in curr_spl.loci.items():
                for method, res in locus.results.items():
                    final_spl.loci[locus_id].results[method] = res
    # Write
    ReportIO.write(final_report, args.output_report)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merge multiple MSI ReportIO from the same samples and same loci. This is used to merge results from multiple classifiers.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-reports', required=True, nargs='+', help='Path to MSI report files (format: JSON).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-report', required=True, help='Path to output report file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    process(args, log)
    log.info("End of job")
