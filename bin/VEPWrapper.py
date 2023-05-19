#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.3'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import logging
import subprocess


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    if len(sys.argv) == 1:
        sys.exit(
            'usage: {} vep_path [vep_arg ...]\n\nWrap the execution of VEP to prevent bug when input file does not contain any variants.\n'.format(
                sys.argv[0]
            )
        )
    vep_cmd = sys.argv[1:]
    ann_tag = vep_cmd[vep_cmd.index("--vcf_info_field") + 1] if "--vcf_info_field" in vep_cmd else "CSQ"
    input_file = vep_cmd[vep_cmd.index("--input_file") + 1]
    output_file = vep_cmd[vep_cmd.index("--output_file") + 1]

    # Logger initialisation
    logging.basicConfig(level=logging.INFO, format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.info("Start with command: {}".format(" ".join(sys.argv)))

    # Check if VCF contains records
    is_empty = True
    with open(input_file) as FH_in:
        for line in FH_in:
            if not line.startswith("#"):
                is_empty = False

    # Process
    if is_empty:  # Create empty output if VCF does not contain any records
        log.info("Skip annotation")
        log.warn('[WARN] The input file "{}" does not contain any records.'.format(input_file))
        with open(input_file) as FH_in:
            with open(output_file, "w") as FH_out:
                for line in FH_in:
                    if line.startswith("#CHROM	POS	ID	REF	ALT	QUAL"):
                        FH_out.write(
                            '##INFO=<ID={},Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature">\n'.format(
                                ann_tag
                            )
                        )
                    FH_out.write(line)
    else:  # Annotate
        log.info("Process annotation")
        subprocess.check_call(vep_cmd)
    log.info("End of job")
