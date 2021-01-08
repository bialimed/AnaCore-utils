#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.annotVcf import AnnotVCFIO, HeaderFormatAttr, HeaderSampleAttr
from anacore.hgvs import HGVSProtChange
from anacore.sv import HashedSVIO
import argparse
import copy
import json
import logging
import os
from owlready2 import get_ontology, owl
import sys


def varCatAreCompatible(variant_cat, evidence_cat):
    return variant_cat.startswith(evidence_cat) or evidence_cat.startswith(variant_cat)


def getVariantCategory(variant):
    variant_category = None
    if "SVTYPE" in record.info:
        cat_by_type = {
            "BND": "fusion",
            "CNV": "mutation",  # "copy_number",
            "DUP:TANDEM": "mutation>insertion",  # "copy_number>gain>tandem",
            "DEL:ME": "copy_number>loss",
            "INS:ME": "copy_number>gain",
            "INV": "inversion",
            "DUP": "mutation>insertion",  # "sv>insertion>duplication",
            "INS": "mutation>insertion",  # "sv>insertion",
            "DEL": "mutation>deletion"  # "sv>deletion"
        }
        variant_category = cat_by_type[record.info["SVTYPE"]]
    else:
        variant_category = variant.type()
        if variant_category == "variation":
            variant_category = "mutation>mnv"
        elif variant_category == "indel":
            variant_category = "mutation>insertion" if variant.isInsertion() else "mutation>deletion"
        else:
            variant_category = "mutation>snv"
    return variant_category


def loadEvidences(in_evidences):
    evidences_by_gene_id = {}
    with HashedSVIO(in_evidences) as reader:
        for record in reader:
            if record["entrez_id"] not in evidences_by_gene_id:
                evidences_by_gene_id[record["entrez_id"]] = list()
            evidences_by_gene_id[record["entrez_id"]].append(record)
    return evidences_by_gene_id


def getDiseaseElt(in_diseases_ontology, disease_id, disease_term):
    if disease_id is not None or disease_term is not None:
        if disease_id is None:
            ontology = get_ontology("file://{}".format(in_diseases_ontology)).load()
            diseases = ontology.search(label=disease_term)
            if len(diseases) == 0:
                raise Exception('The disease "{}" cannot be found in ontology {}.'.format(disease_term, in_diseases_ontology))
            elif len(diseases) != 1:
                raise Exception('The disease "{}" is found several times in ontology {}.'.format(disease_term, in_diseases_ontology))
            disease_id = str(diseases[0].id[0].split(":")[1])
        elif disease_term is None:
            ontology = get_ontology("file://{}".format(in_diseases_ontology)).load()
            disease_term = str(ontology.search(id="DOID:{}".format(disease_id))[0].label[0])
    return disease_id, disease_term


def getDiseaseAncestors(in_diseases_ontology, vcf_handle):
    ontology = get_ontology("file://{}".format(in_diseases_ontology)).load()
    disease_ancestors_by_spl = {}
    for spl_name in vcf_handle.samples:
        spl_doid = vcf_handle.sample_info[spl_name].doid
        if spl_doid is None or spl_doid == "":
            disease_ancestors_by_spl[spl_name] = {}
        else:
            disease = ontology.search(id="DOID:{}".format(spl_doid))[0]
            disease_ancestors_by_spl[spl_name] = {elt.name.replace("DOID_", "") for elt in disease.ancestors() if elt != owl.Thing}
    return disease_ancestors_by_spl


def getEvidences(record, annot_field):
    processed_associations = set()
    retained_evidences = list()
    for annot in record.info[annot_field]:
        if annot["HGVSp"] is not None: ################# and not (annot["HGVSp"].endswith("p.?") or annot["HGVSp"].endswith("=")):  ####################
            annot_association_id = "{}:\t{}".format(annot["SYMBOL"], annot["HGVSp"].split(":", 1)[1])
            if annot_association_id not in processed_associations:  # Skip association on other annotation with the same impact
                processed_associations.add(annot_association_id)
                associations = evidences_by_gene_id[annot["Gene"]] if annot["Gene"] in evidences_by_gene_id else []
                log.debug("{}: {} initial associations".format(record.getName(), len(associations)))
                for curr_asso in associations:
                    if varCatAreCompatible(variant_category, curr_asso["variant_category"]):
                        if curr_asso["HGVSp_change"] == "":
                            retained_evidences.append(curr_asso)
                        else:
                            record_hgvs_p = HGVSProtChange.fromStr(annot["HGVSp"].rsplit("p.", 1)[1])
                            record_hgvs_p.predicted = False
                            if str(record_hgvs_p).startswith(curr_asso["HGVSp_change"]): ###############pb Val15insA et Val15insAV
                                retained_evidences.append(curr_asso)
    return retained_evidences


def addHeaderFields(writer, disease_id=None, disease_term=None):
    writer.format["EVID_PS"] = HeaderFormatAttr("EVID_PS", "The higher clinical evidence level for this precise variant in sample disease.", number=1)
    writer.format["EVID_PA"] = HeaderFormatAttr("EVID_PA", "The higher clinical evidence level for this precise variant in all diseases.", number=1)
    writer.format["EVID_IS"] = HeaderFormatAttr("EVID_IS", "The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in sample disease.", number=1)
    writer.format["EVID_IA"] = HeaderFormatAttr("EVID_IA", "The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in all diseases.", number=1)
    for spl_name in writer.samples:
        if spl_name in writer.sample_info:
            spl_attr = writer.sample_info[spl_name]
            doid = ""
            if "doid" in spl_attr.keys():
                doid = spl_attr.doid
            else:
                spl_attr.case_by_attr["doid"] = "DOID"
            if disease_id:
                doid = disease_id
            spl_attr.doid = doid
            do_term = ""
            if "do_term" in spl_attr.keys():
                do_term = spl_attr.do_term
            else:
                spl_attr.case_by_attr["do_term"] = "DO_term"
            if disease_term:
                do_term = disease_term
            spl_attr.do_term = do_term
        else:
            writer.sample_info[spl_name] = HeaderSampleAttr(
                spl_name,
                doid="" if disease_id is None else disease_id,
                doterm="" if disease_term is None else disease_term
            )
            writer.sample_info[spl_name].case_by_attr["doid"] = "DOID"


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add clinical evidence level to known variants.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    group_disease = parser.add_mutually_exclusive_group()
    group_disease.add_argument('-d', '--disease-id', help='Replace samples disease by this value. It must be a valid ID in disease ontology.')
    group_disease.add_argument('-t', '--disease-term', help='Replace samples disease by this value. It must be a valid TERM in disease ontology.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-e', '--input-evidences', required=True, help='Path to the clinical evidence file from AnaCore-utils/bin/enrichCivic.py (format: TSV).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file. (format: VCF).')
    group_input.add_argument('-n', '--input-diseases-ontology', required=True, help='Path to the diseases ontology file. (format: OWL).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the annotated file. (format: VCF).')
    group_output.add_argument('-s', '--output-evidences', required=True, help='Path to the evidences associated with variants. (format: TSV).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.DEBUG)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_unknown_significance = 0
    nb_with_association = 0
    disease_id, disease_term = getDiseaseElt(args.input_diseases_ontology, args.disease_id, args.disease_term)
    evidences_by_gene_id = loadEvidences(args.input_evidences)
    details_by_spl = {}
    with AnnotVCFIO(args.output_variants, "w") as writer:
        with AnnotVCFIO(args.input_variants) as reader:
            # Header
            writer.copyHeader(reader)
            addHeaderFields(writer, disease_id, disease_term)
            writer.writeHeader()
            disease_ancestors_by_spl = getDiseaseAncestors(args.input_diseases_ontology, writer)
            for spl_name, spl in writer.sample_info.items():
                details_by_spl[spl.id] = {
                    "variants": [],
                    "disease": {"doid": spl.doid, "term": spl.doterm}
                }
            # Records
            for record in reader:
                log.debug(
                    "start {} ({})".format(
                        record.getName(),
                        sorted(set(annot["SYMBOL"] for annot in record.info[args.annotation_field]))
                    )
                )
                # Get corresponding clinical evidences
                variant_category = getVariantCategory(record)
                retained_evidences = getEvidences(record, args.annotation_field)
                # Select best evidences
                best_evidence = {}
                for spl_name in record.samples:
                    best_evidence[spl_name] = {
                        "ps": "",  # precise_variant_same_disease
                        "pa": "",  # precise_variant_all_diseases
                        "is": "",  # imprecise_variant_same_disease
                        "ia": ""  # imprecise_variant_all_diseases
                    }
                    details_by_spl[spl_name]["variants"].append({
                        "gene": ";".join([elt["SYMBOL"] for elt in record.info[args.annotation_field]]),
                        "genomic_alteration": record.getName(),
                        "HGVSp": ";".join([("" if elt["HGVSp"] is None else elt["HGVSp"]) for elt in record.info[args.annotation_field]]),
                        "evidences": []
                    })
                if len(retained_evidences) != 0:
                    nb_with_association += 1
                    for curr_evidence in sorted(retained_evidences, key=lambda elt: elt["evidence_level"]):
                        detailed_annot = {
                            "subject": curr_evidence["cleaned_variant"],
                            "subject_category": curr_evidence["variant_category"],
                            "disease": curr_evidence['disease'],
                            "doid": curr_evidence["doid"],
                            "level": curr_evidence['evidence_level'],
                            "status": curr_evidence['evidence_status'],
                            "type": curr_evidence['evidence_type'],
                            "direction": curr_evidence['evidence_direction'],
                            "drugs": curr_evidence['drugs'],
                            "drug_interaction_type": curr_evidence['drug_interaction_type'],
                            "clinical_significance": curr_evidence['clinical_significance'],
                            "citation": "{}:{}".format(curr_evidence['source_type'], curr_evidence['citation_id']),
                            "source": "CIViC",
                            "specificity": {"variant": "imprecise", "disease": "all"}
                        }
                        for spl_name in record.samples:
                            disease_ancestors = disease_ancestors_by_spl[spl_name]
                            evidence_for_spl = copy.deepcopy(detailed_annot)
                            # Best evidence
                            variant_specificity = ["i"]
                            if curr_evidence["HGVSp_change"] != "":  ########################
                                variant_specificity = ["p", "i"]
                                evidence_for_spl["specificity"]["variant"] = "precise"
                            disease_specificity = ["a"]
                            if curr_evidence['doid'] in disease_ancestors:
                                disease_specificity = ["a", "s"]
                                evidence_for_spl["specificity"]["disease"] = "specific"
                            for curr_v_spe in variant_specificity:
                                for curr_d_spe in disease_specificity:
                                    best_tag = "{}{}".format(curr_v_spe, curr_d_spe)
                                    prev = best_evidence[spl_name][best_tag]
                                    if prev == "" or prev > curr_evidence['evidence_level']:
                                        best_evidence[spl_name][best_tag] = curr_evidence['evidence_level']
                            # Add to detail
                            details_by_spl[spl_name]["variants"][-1]["evidences"].append(evidence_for_spl)
                else:
                    nb_unknown_significance += 1
                # Add best evidence to record
                for spl_name in record.samples:
                    spl_disease = writer.sample_info[spl_name].doid
                    for key, val in best_evidence[spl_name].items():
                        if key[1] == "s" and spl_disease == "":  # Evidence for disease when sample disease is not provided
                            record.samples[spl_name]["EVID_{}".format(key.upper())] = None
                        else:  # Evidence for all diseases or when sample disease is provided
                            record.samples[spl_name]["EVID_{}".format(key.upper())] = val
                record.format.extend(["EVID_PS", "EVID_PA", "EVID_IS", "EVID_IA"])
                # Write record
                writer.write(record)
    #
    # with HashedSVIO(args.output_evidences, "w") as writer_evidences:
    #     writer_evidences.titles = [
    #         "sample", "variant_gene", "variant_hgvsp", "variant_genomic_alteration",
    #         "evidence_variant", "evidence_disease", "evidence_level", "evidence_status",
    #         "evidence_type", "evidence_direction", "drugs", "drug_interaction_type",
    #         "clinical_significance", "citation", "source", "evidence_doid",
    #         "evidence_variant_category", "evidence_variant_specificity",
    #         "evidence_disease_specificity"
    #     ]
    #     for spl_name, spl in details_by_spl.items():
    #         for curr_evidence in spl["evidences"]:
    #             curr_evidence["sample"] = spl_name
    #             writer_evidences.write(curr_evidence)
    with open(args.output_evidences, "w") as writer_evidences:
        json.dump(details_by_spl, writer_evidences)


    log.info("{:.2%} of variants have an association ({}/{})".format(
        nb_with_association / (nb_with_association + nb_unknown_significance),
        nb_with_association,
        nb_with_association + nb_unknown_significance
    ))
    log.info("End of job")
