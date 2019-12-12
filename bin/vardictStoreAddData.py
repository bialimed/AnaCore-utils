#!/usr/bin/env python3

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
from csv import DictReader
from anacore.vcf import VCFIO, HeaderInfoAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def getSupports(variants_tsv):
    """
    Return by variant name the number of reads supporting each allele (ref and alt) on each strand orientation (forward and reverse).

    :param variants_tsv: Path to the variants file outputted by vardict first step (format: TSV).
    :type variants_tsv: str
    :return: By variant name the number of reads supporting each allele (ref and alt) on each strand orientation (forward and reverse).
    :rtype: dict
    """
    supports_by_id = {}
    with open(variants_tsv) as handle:
        titles = [
            "Sample",
            "Gene",
            "Chr",
            "Start",
            "End",
            "Ref",
            "Alt",
            "Depth",
            "AltDepth",
            "RefFwdReads",
            "RefRevReads",
            "AltFwdReads",
            "AltRevReads",
            "Genotype",
            "AF",
            "Bias",
            "PMean",
            "PStd",
            "QMean",
            "QStd",
            "MapQ",
            "QRatio",
            "HiFreq",
            "ExtraAF",
            "Others",  # Change between versions
            # "shift3",
            # "MSI",
            # "MSINT",
            # "NM",
            # "HiCnt",
            # "HiCov",
            # "5pFlankSeq",
            # "3pFlankSeq",
            # "Seg",
            # "VarType"
        ]
        reader = DictReader(handle, delimiter="\t", fieldnames=titles)
        for record in reader:
            id = "{}:{}={}/{}".format(
                record["Chr"],
                record["Start"],
                record["Ref"],
                record["Alt"]
            )
            supports_by_id[id] = {
                "SRF": int(record["RefFwdReads"]),
                "SRR": int(record["RefRevReads"]),
                "SAF": int(record["AltFwdReads"]),
                "SAR": int(record["AltRevReads"])
            }
    return supports_by_id


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add strand support information in vardict VCF from vardict TSV output.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-p', '--input-pre-variants', required=True, help='Path to the variants file outputted by vardict first step (format: TSV).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the vardict variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the outputted variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    supports_by_id = getSupports(args.input_pre_variants)
    with VCFIO(args.output_variants, "w") as FH_out:
        with VCFIO(args.input_variants) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.info["SAF"] = HeaderInfoAttr("SAF", "Number of reads supporting the alternative allele in forward strand.", type="Integer")
            FH_out.info["SAR"] = HeaderInfoAttr("SAR", "Number of reads supporting the alternative allele in reverse strand.", type="Integer")
            FH_out.info["SRF"] = HeaderInfoAttr("SRF", "Number of reads supporting the reference allele in forward strand.", type="Integer")
            FH_out.info["SRR"] = HeaderInfoAttr("SRR", "Number of reads supporting the reference allele in reverse strand.", type="Integer")
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                for key, val in supports_by_id[record.getName()].items():
                    record.info[key] = val
                FH_out.write(record)
    log.info("End of job")
