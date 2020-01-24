#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import sys
import logging
import argparse
from anacore.annotVcf import AnnotVCFIO
from anacore.fusion import FusionFileReader


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Convert TSV output coming from several popular fusions callers to VCF.")
    parser.add_argument('-s', '--sample-name', required=True, help='The sample name.')
    parser.add_argument('-a', '--annotation-field', default="FCANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-fusions', required=True, help='Path to the output of fusion caller (format: TSV).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-fusions', required=True, help='Path to the output file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))
    log.info("Version: " + str(__version__))

    # Process
    with AnnotVCFIO(args.output_fusions, "w", args.annotation_field) as writer:
        with FusionFileReader.factory(args.input_fusions, "r", args.sample_name, args.annotation_field) as reader:
            # Header
            reader.__class__.setVCFHeader(writer, args.annotation_field)
            writer.samples = [args.sample_name]
            writer.writeHeader()
            # Records
            for first_bnd, second_bnd in reader:
                writer.write(first_bnd)
                writer.write(second_bnd)
    log.info("End of job")
