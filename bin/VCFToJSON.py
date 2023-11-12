#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.8.0'

from anacore.annotVcf import AnnotVCFIO, getAlleleRecord
import argparse
import json
import logging
import os
import re
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAnnotSummary(allele_record, initial_alt, annot_field="ANN", pop_prefixes=None, pathogenicity_fields=None, logger=None):
    """
    Return a summary of the diffrent annotations of the variant. This summary is about identical known variants (xref), AF in populations (pop_AF), annotations of the variant and annotations of the collocated variants.

    :param allele_record: The variant.
    :type allele_record: anacore.vcf.VCFRecord
    :param initial_alt: The alternative variant before any transformation (standardization, ...). It must be the alternative variant directly extract from the annotated VCF.
    :type initial_alt: str
    :param annot_field: Field used to store annotations.
    :type annot_field: str
    :param pop_prefixes: The prefixes used to determine database name in population allele frequency fields. Example: the prefix "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF.
    :type pop_prefixes: list
    :param pathogenicity_fields: The titles of fields used to store pathogenicity predictor results. Example: SIFT, PolyPhen, CADD_PHRED.
    :type pathogenicity_fields: str
    :param logger: The logger object.
    :type loggger: logging.Logger
    :return: First the dentical known variants (e.g. {"cosmic": ["COSM14", "COSM15"], "dbSNP":[]}), second AF in populations (e.g. [{"source":"1KG", "name":"Global", "AF":0.85}]), third annotations of the variant and fourth annotations of the collocated variants.
    :rtype: list
    :warnings: The allele_record must only contains one variant.
    """
    xref = {"cosmic": set(), "dbSNP": set(), "HGMD": set(), "Unknown": set()}
    pop_AF = dict()
    variant_annot = list()
    collocated_annot = list()
    for annot in allele_record.info[annot_field]:
        is_self_variant = (initial_alt.upper() == annot["Allele"].upper())
        # Similar knowns variants
        if is_self_variant and annot["Existing_variation"] is not None:
            for db_id in annot["Existing_variation"].split("&"):
                if db_id.startswith("rs"):
                    xref["dbSNP"].add(db_id)
                elif db_id.startswith("COS"):
                    xref["cosmic"].add(db_id)
                elif db_id.startswith("CM"):
                    xref["HGMD"].add(db_id)
                else:
                    xref["Unknown"].add(db_id)
                    if logger is not None:
                        logger.warning('The database using the variant ID "{}" is not managed by "{}".'.format(db_id, sys.argv[0]))
        # Allele frequency in populations
        if is_self_variant:
            for key in annot:
                if key.endswith("_AF") and annot[key] is not None:
                    source, name = getPopInfo(key, pop_prefixes, logger)
                    pop_id = source + "_" + name
                    for curr_AF in annot[key].split("&"):
                        if pop_id not in pop_AF:
                            pop_AF[pop_id] = {
                                "source": source,
                                "name": name,
                                "AF": float(curr_AF)
                            }
                        else:
                            if pop_AF[pop_id]["AF"] != float(curr_AF):
                                raise Exception(
                                    'The allele frequency for the variant {} in population {} is reported several times with different values in "{}".'.format(
                                        allele_record.getName(), pop_id, args.input_variants
                                    )
                                )
        # Annotations
        annot_container = variant_annot if is_self_variant else collocated_annot
        json_annot = {
            "subject": {"symbol": annot["SYMBOL"], "feature": annot["Feature"], "feature_type": annot["Feature_type"]},
            "changes": {"HGVSc": annot["HGVSc"], "HGVSp": annot["HGVSp"]},
            "conseq": annot["Consequence"],
            "pathogenicity": getPathogenicityPredictors(annot, pathogenicity_fields),
            "is_main": annot["PICK"] == "1" if "PICK" in annot else None  # Tag main annotations if flag_pick has been used
        }
        if "EXON" in annot and "INTRON" in annot:  # Add splice_segment only if calculated
            json_annot["pos"] = dict()
            if annot["Feature_type"] == "Transcript":
                if annot["EXON"]:
                    exon_pos = int(annot["EXON"].split("/")[0])
                    json_annot["pos"]["transcript"] = {"exon": exon_pos}
                elif annot["INTRON"]:
                    intron_pos = int(annot["INTRON"].split("/")[0])
                    json_annot["pos"]["transcript"] = {"intron": intron_pos}
                # else: pass  # if variant is up/downstream of transcript (ex: TERT promoter)
        annot_container.append(json_annot)
    xref = {db: list(xref[db]) for db in xref}
    pop_AF = list(pop_AF.values())
    return xref, pop_AF, variant_annot, collocated_annot


def getPathogenicityPredictors(annot, pathogenicity_fields=None):
    """
    Return by predictor the predicted pathogenicity.

    :param annot: The information from an annotation feature.
    :type annot: dict
    :param pathogenicity_fields: The titles of fields used to store pathogenicity predictor results. Example: SIFT, PolyPhen, CADD_PHRED.
    :type pathogenicity_fields: str
    :return: by predictor the predicted pathogenicity.
    :rtype: dict
    """
    pathogenicity_fields = ["CLIN_SIG", "CADD_PHRED", "MetaLR_rankscore", "VEST3_rankscore"] if pathogenicity_fields is None else pathogenicity_fields
    rename = {
        "CLIN_SIG": "ClinVar",
        "CADD_PHRED": "CADD_phred"
    }
    score_by_predictor = {}
    for key in pathogenicity_fields:
        if key in annot and annot[key] != "":
            source = key
            if key in rename:
                source = rename[key]
            score_by_predictor[source] = annot[key]
    return score_by_predictor


def getPopInfo(annot_key, pop_prefixes=None, logger=None):
    """
    Return source and name of a sub-population studied in large genomic programs (1KG, ExAC, ...) from the annotation tag.

    :param annot_key: The tag used in annotation to store AF of the variant in sub-population.
    :type annot_key: str
    :param pop_prefixes: The prefixes used to determine database name in population allele frequency fields (example: "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF).
    :type pop_prefixes: list
    :param logger: The logger object.
    :type loggger: logging.Logger
    :return: The source (project name) and the name of the sub-population.
    :rtype: list
    """
    pop_lower_prefixes = ["exac", "gnomad", "1kg", "esp"] if pop_prefixes is None else [elt.lower() for elt in pop_prefixes]
    source = None
    name = None
    if annot_key == "AF":
        source = "1KG"
        name = "Global"
    elif "_" in annot_key and annot_key.lower().split("_")[0] in pop_lower_prefixes:  # ExAC_AF, ExAC_Adj_AF, ExAC_AFR_AF, ExAC_AMR_AF, ...
        source = annot_key.split("_")[0]
        name = "Global" if annot_key.count('_') == 1 else annot_key.split("_")[1]
    elif annot_key.count('_') == 1:  # AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, ...
        source = "1KG"
        name = annot_key.split("_")[0]
    else:
        if logger is not None:
            logger.warning('The population information stored with tag "{}" cannot be used by "{}".'.format(annot_key, sys.argv[0]))
    return source, name


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Converts VCF annotated with VEP in JSON format.')
    parser.add_argument('-r', '--assembly-id', help='ID of the reference used for the variant calling (example: GRCh38.p12).')
    parser.add_argument('-t', '--pathogenicity-fields', default=["CLIN_SIG", "CADD_PHRED", "MetaLR_rankscore", "VEST3_rankscore"], nargs='+', help='The titles of fields used to store pathogenicity predictor results (example: SIFT, PolyPhen, CADD_PHRED). [Default: %(default)s]')
    parser.add_argument('-p', '--populations-prefixes', default=["gnomadg", "1kg"], nargs='+', help='The prefixes used to determine database name in population allele frequency fields (example: "gnomAD" is used in gnomAD_AF, gnomAD_EUR_AF). [Default: %(default)s]')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-source', default=None, help='Add source of the calling in support information.')
    parser.add_argument('-m', '--merged-sources', action="store_true", help='Indicates that variants file come from a merge of several variants calling with anacore.mergeVCFCaller.py.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the file file containing variants and annotated with VEP v88+ (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='The path to the file outputted file (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Convert VCF in python dict
    json_data = list()
    with AnnotVCFIO(args.input_variants, "r", args.annotation_field) as FH_vcf:
        # Get sources IDs for VCF coming from merged sources
        id_by_src = None
        if args.merged_sources:
            SRC_id_desc = FH_vcf.info["SRC"].description.split("Possible values: ")[1].replace("'", '"')
            id_by_src = json.loads(SRC_id_desc)
        # Records
        for record in FH_vcf:
            for idx_alt, alt in enumerate(record.alt):
                allele_record = getAlleleRecord(FH_vcf, record, idx_alt)
                allele_record.normalizeSingleAllele()
                curr_json = dict()
                # Coord information
                curr_json["coord"] = {
                    "region": allele_record.chrom,
                    "pos": allele_record.pos,
                    "ref": allele_record.ref,
                    "alt": allele_record.alt[0],
                    "assembly": (None if args.assembly_id is None else args.assembly_id)
                }
                # Support information
                if not args.merged_sources:  # The VCF contains results from one variants caller
                    curr_json["supports"] = [{
                        "filters": allele_record.filter,
                        "quality": allele_record.qual,
                        "libraries": [{"alt_depth": allele_record.getAltAD(library)[0], "depth": allele_record.getDP(library), "name": library} for library in FH_vcf.samples],
                        "source": args.calling_source
                    }]
                else:  # The VCF contains results from several variants caller
                    curr_json["supports"] = []
                    for curr_idx, curr_src in enumerate(allele_record.info["SRC"]):
                        src_id = id_by_src[curr_src]
                        src_filters = []
                        for curr_filter in allele_record.filter:
                            if not re.match(r"^s\d+_", curr_filter):  # Filters common between all sources
                                src_filters.append(curr_filter)
                            elif curr_filter.startswith("{}_".format(src_id)):  # Filters coming from the current source
                                src_filters.append(curr_filter.split("_", 1)[1])
                        curr_json["supports"].append({
                            "filters": src_filters,
                            "quality": (allele_record.info["{}_VCQUAL".format(src_id)] if "{}_VCQUAL".format(src_id) in allele_record.info else None),
                            "libraries": [{"alt_depth": library["ADSRC"][curr_idx], "depth": library["DPSRC"][curr_idx], "name": name} for name, library in allele_record.samples.items()],
                            "source": curr_src
                        })
                # Annotations
                if FH_vcf.annot_field in allele_record.info:
                    # Identical known variants, AF in populations and annotations
                    curr_json["xref"], curr_json["pop_AF"], curr_json["annot"], curr_json["collocated_annot"] = getAnnotSummary(
                        allele_record,
                        alt,
                        args.annotation_field,
                        args.populations_prefixes,
                        args.pathogenicity_fields,
                        log
                    )
                json_data.append(curr_json)
                # Evidences
                if "EVID_PA" in FH_vcf.format:
                    curr_json["evidences"] = {}
                    for spl, curr_evidence in allele_record.samples.items():
                        curr_json["evidences"][spl] = {
                            "prec_all_dis": curr_evidence["EVID_PA"],
                            "imp_all_dis": curr_evidence["EVID_IA"],
                        }
                        if curr_evidence["EVID_PS"] is not None:
                            curr_json["evidences"][spl]["prec_same_dis"] = curr_evidence["EVID_PS"]
                            curr_json["evidences"][spl]["imp_same_dis"] = curr_evidence["EVID_IS"]
    # Write output file
    with open(args.output_variants, "w") as FH_out:
        FH_out.write(
            json.dumps(json_data, default=lambda o: o.__dict__, sort_keys=True)
        )
    log.info("End of job")
