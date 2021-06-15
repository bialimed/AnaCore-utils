#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse
from copy import deepcopy
from anacore.vcf import VCFIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def getVariantsByName(callers, haplotyped_variants_files):
    """
    Return variants detection and haplotype by variant name.

    :param callers: List of variants callers corresponding to haplotyped_variants_files.
    :type callers: list
    :param haplotyped_variants_files: Paths to the variants files after merging by haplotype (format: VCF).
    :type haplotyped_variants_files: list
    :return: Variants detection and haplotype by variant name (example: {"chr1:1235448=ATAG/GTAC": {"mutect2": {"merged": True, "sub": {"chr1:1235448=A/G", "chr1:1235451=G/C"}}, "freebayes": {"merged": False, "sub": {"chr1:1235448=ATAG/GTAC"}}}}).
    :rtype: dict
    """
    variants_by_name = {}
    for curr_caller, curr_file in zip(callers, haplotyped_variants_files):
        with VCFIO(curr_file) as FH_in:
            for record in FH_in:
                if record.getName() not in variants_by_name:
                    variants_by_name[record.getName()] = {}
                variants_by_name[record.getName()][curr_caller] = {
                    "merged": "MCO_VAR" in record.info,
                    "sub": set(record.info["MCO_VAR"]) if "MCO_VAR" in record.info else set()
                }
    return variants_by_name


def getVariantsToTransform(callers, variants_by_name):
    """
    Return by caller the uniformised variants(s) replacing the previous variant.

    :param callers: List of variants callers.
    :type callers: list
    :param variants_by_name: Variants detection and haplotype by variant name (example: {"chr1:1235448=ATAG/GTAC": {"mutect2": {"merged": True, "sub": {"chr1:1235448=A/G", "chr1:1235451=G/C"}}, "freebayes": {"merged": False, "sub": {"chr1:1235448=ATAG/GTAC"}}}}).
    :type variants_by_name: dict
    :return: By caller the uniformised variants(s) by previous variant (example: {"mutect2": {"chr1:1235448=A/G": {"chr1:1235448=ATAG/GTAC"}, "chr1:1235451=G/C": {"chr1:1235448=ATAG/GTAC"}}}).
    :rtype: dict
    """
    change_by_caller = {curr_caller: {} for curr_caller in callers}
    for rec_name, info_by_caller in variants_by_name.items():
        if len(info_by_caller) > 1:  # Variant is detected by several callers
            # Determines if variants is composed by a merge of several sub-variants and this composition is different by callers
            contains_merged = False
            cooc_change_merging = False
            sub_variants = None
            for (caller, info) in info_by_caller.items():
                if info["merged"]:
                    contains_merged = True
                if sub_variants is None:  # First caller
                    sub_variants = info["sub"]
                else:
                    if len(info["sub"].symmetric_difference(sub_variants)) != 0:
                        cooc_change_merging = True
            # Add in change the variants composed by a merge of several sub-variants and this composition is different by callers
            if contains_merged and cooc_change_merging:
                # Get retained configuration of the variant
                retained_caller = None
                retained_var = None
                for curr_caller in callers:  # Order determines the configuration priority
                    if retained_var is None and curr_caller in info_by_caller:
                        retained_caller = curr_caller
                        if info_by_caller[curr_caller]["merged"]:
                            retained_var = info_by_caller[curr_caller]["sub"]
                        else:
                            retained_var = {rec_name}
                # Store retained configuration by caller configuration
                for (curr_caller, info) in info_by_caller.items():
                    if curr_caller != retained_caller:
                        init_var = {rec_name} if info["merged"] == 0 else info["sub"]
                        if len(init_var.symmetric_difference(retained_var)) != 0:
                            if info["merged"]:
                                new_var = retained_var - init_var
                                changed_var = init_var - retained_var
                                for sub in changed_var:
                                    change_by_caller[curr_caller][sub] = new_var
                            else:
                                change_by_caller[curr_caller][rec_name] = retained_var
    return change_by_caller


def changeAndWrite(in_vcf, out_vcf, change_by_rec):
    """
    Write uniformised variants.

    :param in_vcf: Path to the variants file before haplotyping (format: VCF).
    :type in_vcf: str
    :param out_vcf: Path to the output variants file.
    :type out_vcf: str
    :param change_by_rec: Uniformised variants(s) by previous variant (example: {"chr1:1235448=A/G": {"chr1:1235448=ATAG/GTAC"}, "chr1:1235451=G/C": {"chr1:1235448=ATAG/GTAC"}}).
    :type change_by_rec: dict
    """
    with VCFIO(out_vcf, "w") as handle_out:
        with VCFIO(in_vcf) as handle_in:
            rec_by_name = {}
            # Manage header
            handle_out.copyHeader(handle_in)
            handle_out.writeHeader()
            # Split/Merge variants
            for record in handle_in:
                if record.getName() not in change_by_rec:
                    if record.getName() not in rec_by_name:
                        rec_by_name[record.getName()] = record
                    elif record.getPopAltAD()[0] > rec_by_name[record.getName()].getPopAltAD()[0]:
                        rec_by_name[record.getName()] = record
                else:  # Must be uniformised
                    if "PGT" in record.format and "PID" in record.format and "PS" in record.format:  # Remove mutect haplotype information
                        for spl_name, spl_info in record.samples.items():
                            del(spl_info["PGT"])
                            del(spl_info["PID"])
                            del(spl_info["PS"])
                        record.format = [elt for elt in record.format if elt not in ["PGT", "PID", "PS"]]
                    # Change variant
                    for curr_retained in change_by_rec[record.getName()]:
                        if curr_retained not in rec_by_name:
                            new_record = deepcopy(record)
                            rec_by_name[curr_retained] = new_record
                            new_record.ref, new_record.alt[0] = curr_retained.split("=")[1].split("/")
                            new_record.pos = int(curr_retained.split("=")[0].split(":")[1])
                        elif record.getPopAltAD()[0] > rec_by_name[curr_retained].getPopAltAD()[0]:
                            new_record = deepcopy(record)
                            rec_by_name[curr_retained] = new_record
                            new_record.ref, new_record.alt[0] = curr_retained.split("=")[1].split("/")
                            new_record.pos = int(curr_retained.split("=")[0].split(":")[1])
            # Write records
            for record in sorted(rec_by_name.values(), key=lambda x: (x.chrom, x.pos, x.refEnd(), x.alt[0])):
                handle_out.write(record)


class LoggerAction(argparse.Action):
    """Manages logger level parameters (The value "INFO" becomes logging.info and so on)."""

    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merge or split co-occuring variants to ensure cohesion between callers.')
    parser.add_argument('-c', '--calling-sources', required=True, nargs='+', help='Name of the source in same order of --inputs-variants. The order of caller determines which configuration between split and merge is retained for variants where callers are not agree.')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--inputs-haplotyped', required=True, nargs='+', help='Paths to the variants files after merging by haplotype (format: VCF). One file by variant caller. The order of caller determines which configuration between split and merge is retained for variants where callers are not agree.')
    group_input.add_argument('-i', '--inputs-variants', required=True, nargs='+', help='Paths to the variants files before merging by haplotype (format: VCF). One file by variant caller. The order of caller determines which configuration between split and merge is retained for variants where callers are not agree.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--outputs-variants', required=True, nargs='+', help='Path to the output variants files (format: VCF).')
    args = parser.parse_args()
    if len(args.inputs_haplotyped) != len(args.inputs_variants) or len(args.inputs_haplotyped) != len(args.calling_sources):
        raise argparse.ArgumentError("{}, {} and {} must contain the same number of elements.".format(
            args.inputs_haplotyped.name,
            args.inputs_variants.name,
            args.calling_sources.name
        ))

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Get variants to change
    change_by_caller = getVariantsToTransform(
        args.calling_sources,
        getVariantsByName(args.calling_sources, args.inputs_haplotyped)
    )

    # Write splitted and merged
    for curr_caller, curr_in, curr_out in zip(args.calling_sources, args.inputs_variants, args.outputs_variants):
        changeAndWrite(curr_in, curr_out, change_by_caller[curr_caller])

    # Log
    for curr_caller in args.calling_sources:
        change_by_rec = change_by_caller[curr_caller] if curr_caller in change_by_caller else {}
        nb_split = 0
        for src_var, new_var in change_by_rec.items():
            log.debug("Change {} to {} from {}.".format(src_var, new_var, curr_caller))
            if len(new_var) > 1:
                nb_split += 1
        log.info("Number of changed variants from {}: {} splitted, {} merged.".format(curr_caller, nb_split, len(change_by_rec) - nb_split))
    log.info("End of job")
