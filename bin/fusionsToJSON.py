#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse
from copy import deepcopy
from anacore.fusion import BreakendVCFIO, getStrand


########################################################################
#
# FUNCTIONS
#
########################################################################
def getBreakendInfo(record, annot_field="ANN", assembly_id=None):
    coordinates = {
        "region": record.chrom,
        "pos": record.pos,
        "annot_pos": record.info["ANNOT_POS"],
        "ref": record.ref,
        "alt": record.alt[0],
        "assembly": (None if assembly_id is None else assembly_id),
        "strand": getStrand(record)
    }
    features = {}
    for idx, feature in enumerate(record.info[annot_field]):
        cp_feature = deepcopy(feature)
        del(cp_feature["IN_FRAME"])
        features[str(idx)] = cp_feature
    return {"coordinates": coordinates, "features_by_id": features}


def getSupportMergedSources(record, id_by_src):
    support = {}
    for curr_idx, curr_src in enumerate(record.info["SRC"]):
        src_id = id_by_src[curr_src]
        support[curr_src] = {
            "quality": (record.info["{}_VCQUAL".format(src_id)] if "{}_VCQUAL".format(src_id) in record.info else None),
            "libraries": [
                {
                    "spanning_reads": sample["SRSRC"][curr_idx],
                    "spanning_pairs": sample["PRSRC"][curr_idx],
                    "name": name
                } for name, sample in record.samples.items()
            ]
        }
    return support


def getSupportUniqSource(record, calling_source="unknown"):
    return {
        calling_source: {
            "quality": record.qual,
            "libraries": [
                {
                    "spanning_reads": sample["SR"],
                    "spanning_pairs": sample["PR"],
                    "name": name
                } for name, sample in record.samples.items()
            ]
        }
    }


def getKnownPartners(record, known_field="known_partners"):
    db_by_partners = {}
    for partners_info in record.info["known_partners"]:  # partners_info example: BCR_@_ABL1=cosmic_91:1743,1745|chimerdb_pub-V4:3427,3428
        partners, databases = partners_info.split("=", 1)
        databases = databases.split("|")
        partners_res = {}
        for curr_db in databases:
            db_name, db_entries = curr_db.split(":", 1)
            partners_res[db_name] = db_entries.split(",")
        db_by_partners[partners] = partners_res
    return db_by_partners


def getFusionAnnot(record, mate, annot_field="ANN"):
    fusion_annot = list()
    if len(record.info[annot_field]) == 0:
        if len(mate.info[annot_field]) != 0:  # record: no annot ; mate: annot
            for mate_idx, mate_annot in enumerate(mate.info[annot_field]):
                fusion_annot.append({
                    "features_annotation": [None, str(mate_idx)],
                    "inframe": "."
                })
    else:
        if len(mate.info[annot_field]) == 0:  # record: annot ; mate: no annot
            for record_idx, record_annot in enumerate(record.info[annot_field]):
                fusion_annot.append({
                    "features_annotation": [str(record_idx), None],
                    "inframe": "."
                })
        else:  # record: annot ; mate: annot
            for record_idx, record_annot in enumerate(record.info[annot_field]):
                inframe_by_partner = {}
                for elt in record_annot["IN_FRAME"].split("&"):
                    partner, inframe = elt.split(":")
                    if inframe == "1":
                        inframe = True
                    elif inframe == "0":
                        inframe = False
                    inframe_by_partner[partner] = inframe
                for mate_idx, mate_annot in enumerate(mate.info[annot_field]):
                    fusion_annot.append({
                        "features_annotation": [str(record_idx), str(mate_idx)],
                        "inframe": inframe_by_partner[mate_annot["Feature"]]
                    })
    return fusion_annot


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Converts fusions VCF annotated with annotBND in JSON format.')
    parser.add_argument('-r', '--assembly-id', help='ID of the reference used for the variant calling (example: GRCh38.p12).')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-source', default=None, help='Add source of the calling in support information.')
    parser.add_argument('-m', '--merged-sources', action="store_true", help='Indicates that variants file come from a merge of several fusions callers.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file file containing fusions annotated by AnaCore-utils/annotBND.py (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='The path to the file outputted file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Convert VCF to python dict
    json_data = list()
    with BreakendVCFIO(args.input_variants, "r", args.annotation_field) as reader:
        # Get sources IDs for VCF coming from merged sources
        id_by_src = None
        if args.merged_sources:
            SRC_id_desc = reader.info["SRC"].description.split("Possible values: ")[1].replace("'", '"')
            id_by_src = json.loads(SRC_id_desc)
        # Records
        for record, mate in reader:
            curr_json = dict()
            # Coord information
            curr_json["breakends"] = [
                getBreakendInfo(record, args.annotation_field, args.assembly_id),
                getBreakendInfo(mate, args.annotation_field, args.assembly_id)
            ]
            # Filters
            curr_json["filters"] = record.filter
            # Support information
            if not args.merged_sources:  # The VCF contains results from one variants caller
                ref_source = args.calling_source if args.calling_source is not None else "unknown"
                curr_json["supports"] = {
                    "selected": ref_source,
                    "by_src": getSupportUniqSource(record, ref_source)
                }
            else:  # The VCF contains results from several variants caller
                curr_json["supports"] = {
                    "selected": record.info["REFSRC"],
                    "by_src": getSupportMergedSources(record, id_by_src)
                }
            # Fusion level annotations
            if reader.annot_field in record.info or reader.annot_field in mate.info:
                curr_json["annotations"] = getFusionAnnot(record, mate, args.annotation_field)
            # Known partners
            if "known_partners" in record.info:
                curr_json["known_partners"] = getKnownPartners(record, "known_partners")
            json_data.append(curr_json)

    # Write output file
    with open(args.output_variants, "w") as writer:
        writer.write(
            json.dumps(json_data, default=lambda o: o.__dict__, sort_keys=True)
        )
    log.info("End of job")
