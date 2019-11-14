#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import pysam
import argparse
from anacore.illumina import getInfFromSeqID


########################################################################
#
# FUNCTIONS
#
########################################################################
def getGatkRgId(read_id):
    """
    Return the read group ID of the read. GATK RG ID format: flowcell_id.lane_id.

    :param read_id: The read ID.
    :type read_id: str
    :return: Read group ID
    :rtype: str
    """
    read_info = getInfFromSeqID(read_id)
    return "{}.{}".format(read_info["flowcell_id"], read_info["lane_id"])


def getRGInfo(input_aln, force_spl=None, force_lib=None, barcode=None):
    """
    Return new RG info by uid ([old_rg_id.]flowcell_id.lane_id[.barcode]) and old read RG info by old RG_ID.

    :param input_aln: Path to the previous alignment file (format: BAM).
    :type input_aln: str
    :param force_spl: Force sample name in new RG. Otherwise the samples names are retrieved from the previous reads groups if they exist.
    :type force_spl: str
    :param force_lib: Force library name in new RG. Otherwise the libraries names are retrieved from the previous reads groups if they exist.
    :type force_lib: str
    :param barcode: barcode added in new platform unit.
    :type barcode: str
    :return: new RG info by uid ([old_rg_id.]flowcell_id.lane_id[.barcode]) and old read RG info by old RG_ID.
    :rtype: (dict, dict)
    """
    RG_by_uid = {}
    old_RG_by_rgid = {}
    next_id = 1
    with pysam.AlignmentFile(input_aln, "rb") as FH_in:
        # Header
        header = FH_in.header
        if "RG" in header:
            for group in header["RG"]:
                old_RG_by_rgid[group["ID"]] = group
        # Records
        for curr_read in FH_in.fetch(until_eof=True):
            old_rg = {} if not curr_read.has_tag("RG") else old_RG_by_rgid[curr_read.get_tag("RG")]
            id = getGatkRgId(curr_read.query_name)
            pu = id if barcode is None else "{}.{}".format(id, barcode)
            uid = pu if "ID" not in old_rg else "{}.{}".format(old_rg["ID"], pu)
            if uid not in RG_by_uid:
                RG_by_uid[uid] = {
                    "ID": str(next_id),
                    "PL": "ILLUMINA",
                    "PU": pu
                }
                next_id += 1
                # SM
                if force_spl is not None:
                    RG_by_uid[uid]["SM"] = force_spl
                elif "SM" in old_rg:
                    RG_by_uid[uid]["SM"] = old_rg["SM"]
                # LB
                if force_lib is not None:
                    RG_by_uid[uid]["LB"] = force_lib
                elif "LB" in old_rg:
                    RG_by_uid[uid]["LB"] = old_rg["LB"]
    return RG_by_uid, old_RG_by_rgid


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add RG on reads. ID, PL and PU are determined from reads and optional option barcode ; LB, SM can be retrieved from previous BAM or forced.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_RG = parser.add_argument_group('Reads groups')  # Reads groups
    group_RG.add_argument('-l', '--library-name', help='Force library name. Otherwise the libraries names are retrieved from the previous reads groups if they exist.')
    group_RG.add_argument('-b', '--barcode', help='Add barcode in platform unit.')
    group_RG.add_argument('-s', '--sample-name', help='Force sample name. Otherwise the samples names are retrieved from the previous reads groups if they exist.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='The path to the alignments file (format: BAM).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-aln', required=True, help='The path to the outputted alignments file (format: BAM).')
    args = parser.parse_args()

    # Get the reads groups
    RG_by_uid, old_RG_by_rgid = getRGInfo(args.input_aln, args.sample_name, args.library_name, args.barcode)

    # Change read groups
    with pysam.AlignmentFile(args.input_aln, "rb") as FH_in:
        # Replace RG in header
        new_header = FH_in.header.to_dict()
        new_header["RG"] = list(RG_by_uid.values())
        # Replace RG in reads
        with pysam.AlignmentFile(args.output_aln, "wb", header=new_header) as FH_out:
            for curr_read in FH_in.fetch(until_eof=True):
                old_rg = {} if not curr_read.has_tag("RG") else old_RG_by_rgid[curr_read.get_tag("RG")]
                id = getGatkRgId(curr_read.query_name)
                pu = id if args.barcode is None else "{}.{}".format(id, args.barcode)
                uid = pu if "ID" not in old_rg else "{}.{}".format(old_rg["ID"], pu)
                curr_read.set_tag("RG", RG_by_uid[uid]["ID"])
                FH_out.write(curr_read)
