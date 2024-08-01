#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'

from anacore.vcf import VCFIO, HeaderInfoAttr
import argparse
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def getCleanningRules(variant_caller):
    """
    Return by INFO tag the correct declaration for header and the function to clean the values of this tag in records.

    :param variant_caller: The variant caller used to produce the VCF to fix.
    :type variant_caller: str
    :return: By INFO tag the correct declaration for header and the function to clean the values of this tag in records.
    :rtype: dict
    """
    info_by_caller = {
        "vardict": {
            "REFBIAS": {
                "declaration": HeaderInfoAttr("REFBIAS", "Reference depth by strand", type="Integer", number="2"),
                "process": lambda val: [int(elt) for elt in val.split(":")]
            },
            "VARBIAS": {
                "declaration": HeaderInfoAttr("VARBIAS", "Variant depth by strand", type="Integer", number="2"),
                "process": lambda val: [int(elt) for elt in val.split(":")]
            }
        }
    }
    return info_by_caller[variant_caller]


def getDelInfo(variant_caller):
    """
    Return tags to delete in INFO field.

    :param variant_caller: The variant caller used to produce the VCF to fix.
    :type variant_caller: str
    :return: Tags to delete in INFO field.
    :rtype: set
    """
    del_by_caller = {
        "vardict": ["SAMPLE"]
    }
    return del_by_caller[variant_caller]


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Fix bug in INFO fields format of variant caller outputs.')
    parser.add_argument('-c', '--variant-caller', default="vardict", choices=["vardict"], help='The variant caller used to produce the VCF to fix. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
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
    clean_info = getCleanningRules(args.variant_caller)
    del_info = getDelInfo(args.variant_caller)
    with VCFIO(args.output_variants, "w") as writer:
        with VCFIO(args.input_variants) as reader:
            # Header
            writer.copyHeader(reader)
            for tag in writer.info:
                if tag in clean_info:
                    prev = writer.info[tag]
                    new = clean_info[tag]["declaration"]
                    if prev.type != new.type or prev.number != new.number:
                        writer.info[tag] = new
                    else:
                        del(clean_info[tag])
            for tag in del_info:
                del(writer.info[tag])
            writer.writeHeader()
            # Records
            for record in reader:
                for tag, value in record.info.items():
                    if tag in clean_info:
                        record.info[tag] = clean_info[tag]["process"](value)
                writer.write(record)
    log.info("End of job")
