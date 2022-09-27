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
from anacore.bed import BEDIO
from anacore.sv import HashedSVIO
from anacore.msi.base import Status
from anacore.msi.annot import MSIAnnot


########################################################################
#
# FUNCTIONS
#
########################################################################
def process(args, log):
    """
    Convert MSI status file (splA<tab>status_locus_1<tab>status_locus_2) in MSI annotation file.

    :param args: The namespace extracted from the script arguments.
    :type args: Namespace
    :param log: The logger of the script.
    :type log: loggin.Logger
    """
    # Get targeted loci IDs and names
    loci_in_bed = []
    id_by_name = {}
    with BEDIO(args.input_targets) as FH_in:
        for record in FH_in:
            id = "{}:{}-{}".format(record.chrom, record.start - 1, record.end)
            id_by_name[record.name] = id
            loci_in_bed.append(id)
    if not args.locus_id:
        loci_in_bed = sorted(id_by_name.keys())
    # Write annotation file
    with HashedSVIO(args.input_status, title_starter="") as FH_in:
        loci_in_status = set([elt for elt in FH_in.titles if elt != "sample"])
        if len(set(loci_in_bed) - loci_in_status) > 0:
            msg = "The following loci are defined in targets but are missing from status file: {}".format(set(loci_in_bed) - loci_in_status)
            log.error(msg)
            raise Exception(msg)
        with MSIAnnot(args.output_annotations, "w") as FH_out:
            for record in FH_in:
                for locus in loci_in_bed:
                    if record[locus] not in Status.authorizedValues():
                        msg = 'The status "{}" of the locus {} in sample {} is invalid. It must be: {}'.format(
                            record[locus],
                            locus,
                            record["sample"],
                            Status.authorizedValues()
                        )
                        log.error(msg)
                        raise Exception(msg)
                    FH_out.write({
                        "sample": record["sample"],
                        "locus_position": locus if args.locus_id else id_by_name[locus],
                        "method_id": "model",
                        "key": "status",
                        "value": record[locus],
                        "type": "str"
                    })


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Convert MSI status file (splA<tab>status_locus_1<tab>status_locus_2) in MSI annotation file.')
    parser.add_argument('-l', '--locus-id', action='store_true', help='Use this option if the titles in input-status correspond to the loci IDs (format: chr:start-end with start 0-based) and not the names defined in input-target.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-t', '--input-targets', required=True, help='Path to the file containing locations of microsatellites (format: BED).')
    group_input.add_argument('-s', '--input-status', required=True, help='Path to the file containing for each sample the status of each locus (format: TSV). The title line must contain "sample<tab>loci_1_name...<tab>loci_n_name". Each row defined one sample with the format: "spl_name<tab>loci_1_status<tab>...<tab>loci_n_status<tab>". The status must be "MSS", "MSI" or "Undetermined".')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-a', '--output-annotations', required=True, help='The path to the output file (format: MSIAnnot).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    process(args, log)
    log.info("End of job")
