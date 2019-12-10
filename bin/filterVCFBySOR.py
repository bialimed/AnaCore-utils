#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import math
import logging
import argparse
from anacore.vcf import VCFIO, HeaderFilterAttr, HeaderInfoAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def strandOddRatio(ref_fwd, ref_rev, alt_fwd, alt_rev):
    """
    Return the strand symmetric odds ratio.

    Calculation method come from GATK (details on: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_StrandOddsRatio.php)
    :param ref_fwd: Number of reads supporting the reference allele in forward strand.
    :type ref_fwd: int
    :param ref_rev: Number of reads supporting the reference allele in reverse strand.
    :type ref_rev: int
    :param alt_fwd: Number of reads supporting the alternative allele in forward strand.
    :type alt_fwd: int
    :param alt_rev: Number of reads supporting the alternative allele in reverse strand.
    :type alt_rev: int
    :return: Strand symmetric odds ratio.
    :rtype: float
    """
    r = ((ref_fwd + 1) * (alt_rev + 1)) / ((alt_fwd + 1) * (ref_rev + 1))
    symetrical_ratio = r + (1 / r)
    ref_ratio = (min(ref_fwd, ref_rev) + 1) / (max(ref_fwd, ref_rev) + 1)
    alt_ratio = (min(alt_fwd, alt_rev) + 1) / (max(alt_fwd, alt_rev) + 1)
    strand_odd_ratio = math.log(symetrical_ratio) + math.log(ref_ratio) - math.log(alt_ratio)
    return strand_odd_ratio


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter variants on strand bias estimated by the symmetric odds ratio. Calculation details on https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_annotator_StrandOddsRatio.php.')
    parser.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]')
    group_filter = parser.add_argument_group('Filter')  # Filter
    group_filter.add_argument('-b', '--bias-tag', default="strandRatioBias", help='Tag added to filter field in tag mode. [Default: %(default)s]')
    group_filter.add_argument('-im', '--indel-max-SOR', default=10, type=float, help='Maximum strand symmetric odds ratio for in/del variants. [Default: %(default)s]')
    group_filter.add_argument('-sm', '--substit-max-SOR', default=3, type=float, help='Maximum strand symmetric odds ratio for substitution(s) variants. [Default: %(default)s]')
    group_calculation = parser.add_argument_group('Calculation')  # Calculation
    group_calculation.add_argument('-s', '--SOR-tag', default="SOR", help='Key of the field that will use to store the strand symmetric odds ratio. [Default: %(default)s]')
    group_calculation.add_argument('-rf', '--ref_fwd_tag', default="SRF", help='Key of the field containing the number of reads supporting the reference allele in forward strand. [Default: %(default)s]')
    group_calculation.add_argument('-rr', '--ref_rev_tag', default="SRR", help='Key of the field containing the number of reads supporting the reference allele in reverse strand. [Default: %(default)s]')
    group_calculation.add_argument('-af', '--alt_fwd_tag', default="SAF", help='Key of the field containing the number of reads supporting the alternative allele in forward strand. [Default: %(default)s]')
    group_calculation.add_argument('-ar', '--alt_rev_tag', default="SAR", help='Key of the field containing the number of reads supporting the alternative allele in reverse strand. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', help='Path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', help='Path to the file outputted file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_variants = 0
    nb_filtered = 0
    with VCFIO(args.input_variants) as FH_in:
        with VCFIO(args.output_variants, "w") as FH_out:
            # Header
            FH_out.copyHeader(FH_in)
            FH_out.info[args.SOR_tag] = HeaderInfoAttr(args.SOR_tag, "Strand bias estimated by the symmetric odds ratio test.", type="Float")
            FH_out.filter[args.bias_tag] = HeaderFilterAttr(args.bias_tag, "Strand ratio bias (estimated by the symmetric odds ratio test): substit SOR > {}, InDel SOR > {}.".format(args.substit_max_SOR, args.indel_max_SOR))
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                if len(record.alt) > 1:
                    raise Exception("The multi-allelic variants cannot be processed: {}.".format(record.getName()))
                nb_variants += 1
                is_filtered = False
                # Compute SOR
                record.info[args.SOR_tag] = strandOddRatio(
                    record.info[args.ref_fwd_tag],
                    record.info[args.ref_rev_tag],
                    record.info[args.alt_fwd_tag],
                    record.info[args.alt_rev_tag]
                )
                # Evaluate filter
                if record.type() == "indel":  # InDel
                    if record.info[args.SOR_tag] > args.indel_max_SOR:
                        is_filtered = True
                        nb_filtered += 1
                elif record.info[args.SOR_tag] > args.substit_max_SOR:  # Substitution
                    is_filtered = True
                    nb_filtered += 1
                # Filter record
                if args.mode == "remove":
                    if not is_filtered:
                        if record.filter is None or len(record.filter) == 0:
                            record.filter = ["PASS"]
                        FH_out.write(record)
                else:
                    filters = set()
                    if record.filter is not None and len(record.filter) != 0 and record.filter[0] != "PASS":
                        filters = set(record.filter)
                    if is_filtered:
                        filters.add(args.bias_tag)
                    record.filter = sorted(filters)
                    if len(record.filter) == 0:
                        record.filter = ["PASS"]
                    FH_out.write(record)

    # Log process
    log.info(
        "{:.2%} of variants have been {} ({}/{})".format(
            0 if nb_variants == 0 else nb_filtered / nb_variants,
            "tagged" if args.mode == "tag" else "removed",
            nb_filtered,
            nb_variants
        )
    )
    log.info("End of job")
