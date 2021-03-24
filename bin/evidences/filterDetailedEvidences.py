#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sv import HashedSVIO
from anacore.vcf import VCFIO
import argparse
import json
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def getKeptVariants(in_variants):
    """
    Return the set of variants names by sample from a variants file.

    :param in_variants: Path to the variants file (format: VCF).
    :type in_variants: str
    :return: The set of variants names by sample.
    :rtype: dict
    """
    kept_variants_by_spl = {}
    with VCFIO(in_variants) as reader:
        for spl in reader.samples:
            kept_variants_by_spl[spl] = set()
        for record in reader:
            kept_variants_by_spl[spl].add(record.getName())
    return kept_variants_by_spl


def filterJSON(in_evidences, out_evidences, kept_variants_by_spl):
    """
    Write the filtered version of the detailed evidences file.

    :param in_evidences: Path to the detailed evidences file produced by annotEvidences.py (format: JSON).
    :type in_evidences: str
    :param out_evidences: Path to the filtered evidences file (format: JSON).
    :type out_evidences: str
    :param kept_variants_by_spl: The set of variants names for kept variants by sample.
    :type kept_variants_by_spl: dict
    """
    # Load
    evidences_by_spl = None
    with open(in_evidences) as reader:
        evidences_by_spl = json.load(reader)
    # Filter
    evidence_samples = list(evidences_by_spl.keys())
    for spl in evidence_samples:
        if spl not in kept_variants_by_spl:
            del(evidences_by_spl[spl])
        else:
            variants = evidences_by_spl[spl]["variants"]
            nb_variants = len(variants)
            for variant_idx in reversed(range(nb_variants)):
                variant_name = variants[variant_idx]["genomic_alteration"]
                if variant_name not in kept_variants_by_spl[spl]:
                    del(variants[variant_idx])
    # Write
    with open(out_evidences, "w") as writer:
        json.dump(evidences_by_spl, writer)


def filterTSV(in_evidences, out_evidences, kept_variants_by_spl):
    """
    Write the filtered version of the detailed evidences file.

    :param in_evidences: Path to the detailed evidences file produced by annotEvidences.py (format: TSV).
    :type in_evidences: str
    :param out_evidences: Path to the filtered evidences file (format: TSV).
    :type out_evidences: str
    :param kept_variants_by_spl: The set of variants names for kept variants by sample.
    :type kept_variants_by_spl: dict
    """
    with HashedSVIO(in_evidences) as reader:
        with HashedSVIO(out_evidences, "w") as writer:
            writer.titles = reader.titles
            writer.writeHeader()
            for record in reader:
                spl = record["sample"]
                variant = record["variant_genomic_alteration"]
                if spl in kept_variants_by_spl and variant in kept_variants_by_spl[spl]:
                    writer.write(record)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter detailed evidences produced by annotEvidences.py to correspond to the filtered VCF.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-evidences', required=True, help='Path to the detailed evidences file produced by annotEvidences.py (format: JSON or TSV).')
    group_input.add_argument('-a', '--input-variants', required=True, help='Path to variants file annotated by annotEvidences.py and post-filtered (format: VCF).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-evidences', help='Path to the filtered evidences file (format: TSV or JSON depends on input).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    kept_variants_by_spl = getKeptVariants(args.input_variants)
    is_json = False
    with open(args.input_evidences) as reader:
        is_json = reader.read(1) == "{"
    if is_json:
        filterJSON(args.input_evidences, args.output_evidences, kept_variants_by_spl)
    else:
        filterTSV(args.input_evidences, args.output_evidences, kept_variants_by_spl)
