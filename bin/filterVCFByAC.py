#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import argparse
from anacore.vcf import VCFIO, HeaderFilterAttr


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Filter variants on their AD, AF and DP.')
    parser.add_argument('-m', '--mode', default="tag", choices=["tag", "remove"], help='Select the filter mode. In mode "tag" if the variant does not fit criteria a tag "lowAD" and/or "lowAF" is added in FILTER field. In mode "remove" if the variant does not fit criteria it is removed from the output. [Default: %(default)s]')
    parser.add_argument('-d', '--min-AD', default=4, type=int, help='Filter variants with AD <= than this values. [Default: %(default)s]')
    parser.add_argument('-f', '--min-AF', default=0.02, type=float, help='Filter variants with AF <= than this values. [Default: %(default)s]')
    parser.add_argument('-p', '--min-DP', default=20, type=int, help='Filter variants with DP <= than this values. [Default: %(default)s]')
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
            FH_out.filter["lowAF"] = HeaderFilterAttr("lowAF", "Variants with population AF <= {}".format(args.min_AF))
            FH_out.filter["lowAD"] = HeaderFilterAttr("lowAD", "Variants with population AD <= {}".format(args.min_AD))
            FH_out.filter["lowDP"] = HeaderFilterAttr("lowDP", "Variants with population DP <= {}".format(args.min_DP))
            if "ADSRC" in FH_in.format:  # Variants file come from merging of different calling
                for curr in ["lowAF", "lowAD", "lowDP"]:
                    FH_out.filter[curr].description = FH_out.filter[curr].description + " in all variant calling sources"
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                if len(record.alt) > 1:
                    raise Exception("The multi-allelic variants cannot be processed: {}.".format(record.getName()))
                nb_variants += 1
                new_filters = set()
                if "ADSRC" in FH_in.format:  # Variants file come from merging of different calling
                    pop_AD = [0 for src in record.info["SRC"]]  # Population AD for each calling source
                    pop_DP = [0 for src in record.info["SRC"]]  # Population DP for each calling source
                    for spl_name, spl_info in record.samples.items():
                        for idx, (AD, DP) in enumerate(zip(spl_info["ADSRC"], spl_info["DPSRC"])):
                            pop_AD[idx] += AD
                            pop_DP[idx] += DP
                    pop_AF = [AD / DP for AD, DP in zip(spl_info["ADSRC"], spl_info["DPSRC"])]
                    if max(pop_AF) < args.min_AF:
                        new_filters.add("lowAF")
                    if max(pop_AD) < args.min_AD:
                        new_filters.add("lowAD")
                    if max(pop_DP) < args.min_DP:
                        new_filters.add("lowDP")
                else:
                    if record.getPopAltAF()[0] < args.min_AF:
                        new_filters.add("lowAF")
                    if record.getPopAltAD()[0] < args.min_AD:
                        new_filters.add("lowAD")
                    if record.getPopDP() < args.min_DP:
                        new_filters.add("lowDP")
                if len(new_filters) > 0:
                    nb_filtered += 1
                # Filter record
                if args.mode == "remove":
                    if len(new_filters) == 0:
                        if record.filter is None or len(record.filter) == 0:
                            record.filter = ["PASS"]
                        FH_out.write(record)
                else:
                    old_filters = set()
                    if record.filter is not None and len(record.filter) != 0 and record.filter[0] != "PASS":
                        old_filters = set(record.filter)
                    record.filter = sorted(old_filters | new_filters)
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
