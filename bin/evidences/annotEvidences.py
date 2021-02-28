#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.annotVcf import AnnotVCFIO, HeaderFormatAttr, HeaderSampleAttr
from anacore.hgvs import HGVSProtChange
from anacore.region import Region, RegionList
from anacore.sv import HashedSVIO
import argparse
import copy
import json
import logging
import os
from owlready2 import get_ontology, owl
import sys


def varCatAreCompatible(variant_cat, evidence_cat):
    """
    Return True if the evidence feature category correspond to the variant category.

    :param variant_cat: Category of variant (see getVariantCategory()).
    :type variant_cat: str
    :param evidence_cat: Category of evidence feature.
    :type evidence_cat: str
    :return: True if the evidence feature category correspond to the variant category.
    :rtype: bool
    """
    return variant_cat.startswith(evidence_cat) or evidence_cat.startswith(variant_cat)


def getVariantCategory(variant):
    """
    Return category of variant:
      copy_number or copy_number>gain or copy_number>loss or
      expression or expression>overexpression or expression>underexpression or
      fusion or
      inversion or
      methylation or methylation>demethylation or methylation>hypermethylation or
      mutation or mutation>deletion or mutation>insertion or mutation>substitution or
      phosphorylation

    :param variant: The variant.
    :type variant: anacore.vcf.VCFRecord
    :return: Category of variant.
    :rtype: str
    """
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
            variant_category = "mutation>substitution"
        elif variant_category == "indel":
            variant_category = "mutation>insertion" if variant.isInsertion() else "mutation>deletion"
        else:
            variant_category = "mutation>substitution"
    return variant_category


def loadEvidences(in_evidences):
    """
    Return the list of evidences by gene ENTREZ id.

    :param in_evidences: Path to the clinical evidences file from AnaCore-utils/bin/enrichCivic.py (format: TSV).
    :type in_evidences: str
    :retrun: List of evidences by gene ENTREZ id.
    :rtype: dict
    """
    evidences_by_gene_id = {}
    with HashedSVIO(in_evidences) as reader:
        for record in reader:
            if record["entrez_id"] not in evidences_by_gene_id:
                evidences_by_gene_id[record["entrez_id"]] = list()
            evidences_by_gene_id[record["entrez_id"]].append(record)
    return evidences_by_gene_id


def getOntologyElementByID(ontology, disease_id):
    """
    Return disease element from ontology.

    :param ontology: Disease ontology.
    :type ontology: owlready2.namespace.Ontology
    :param disease_id: The disease ID in ontology.
    :type disease_id: str
    :return: The disease element from ontology.
    :rtype: owlready2.entity.ThingClass
    """
    element = None
    diseases = ontology.search(id=disease_id)
    if len(diseases) == 0:
        diseases = ontology.search(hasAlternativeId=disease_id)
        if len(diseases) > 1:
            raise Exception('The disease with ID "{}" is found several times in ontology.'.format(disease_id))
    elif len(diseases) > 1:
        raise Exception('The disease with ID "{}" is found several times in ontology.'.format(disease_id))
    if len(diseases) == 1:
        element = diseases[0]
    return element


def getOntologyElementByTerm(ontology, disease_term, case_sensitive=True):
    """
    Return disease element from ontology.

    :param ontology: Disease ontology.
    :type ontology: owlready2.namespace.Ontology
    :param disease_term: The disease term or exact synonym in ontology.
    :type disease_term: str
    :param case_sensitive: The search of disease term is case sensitive.
    :type case_sensitive: bool
    :return: The disease element from ontology.
    :rtype: owlready2.entity.ThingClass
    """
    element = None
    diseases = ontology.search(label=disease_term, _case_sensitive=case_sensitive)
    if len(diseases) == 0:
        diseases = ontology.search(hasExactSynonym=disease_term, _case_sensitive=case_sensitive)
        if len(diseases) > 1:
            raise Exception('The disease "{}" is found several times in ontology.'.format(disease_term))
    elif len(diseases) > 1:
        raise Exception('The disease "{}" is found several times in ontology.'.format(disease_term))
    if len(diseases) == 1:
        element = diseases[0]
    return element


def getDiseaseElt(ontology, disease_id, disease_term):
    """
    Return disease id and term from ontology if one of this element is previously defined.

    :param ontology: Disease ontology.
    :type ontology: owlready2.namespace.Ontology
    :param disease_id: The disease ID in ontology.
    :type disease_id: str
    :param disease_term: The disease term or exact synonym in ontology.
    :type disease_term: str
    :return: The tuple disease id and disease term from ontology.
    :rtype: (str, str)
    """
    if disease_id is not None or disease_term is not None:
        if disease_id is None:
            disease = getOntologyElementByTerm(ontology, disease_term, False)
            if disease is None:
                raise Exception('The disease "{}" cannot be found in ontology.'.format(disease_term))
            disease_id = str(disease.id[0])
            disease_term = str(disease.label[0])
        elif disease_term is None:
            disease = getOntologyElementByID(ontology, disease_id)
            if disease is None:
                raise Exception('The disease with DOID "{}" cannot be found in ontology.'.format(disease_id))
            disease_term = str(disease.label[0])
            disease_id = str(disease.id[0])
    return disease_id, disease_term


def getDiseaseAncestors(ontology, vcf_handle):
    """
    Return by sample the list containing his disease and the more generic terms including the disease. Example of list for sample with gastric cardia carcinoma: ["gastric cardia carcinoma", "stomach carcinoma", "cardia cancer", "stomach cancer", "stomach disease", "gastrointestinal system cancer", "gastrointestinal system disease", "cancer"]).

    :param ontology: Disease ontology.
    :type ontology: owlready2.namespace.Ontology
    :vcf_handle: Handle to the VCF containing samples markup with DOID attribute in header.
    :type vcf_handle: anacore.vcf.VCFIO
    :return: by sample the list of his disease and their disesases parents.
    :rtype: dict
    """
    disease_ancestors_by_spl = {}
    for spl_name in vcf_handle.samples:
        spl_doid = vcf_handle.sample_info[spl_name].doid
        if spl_doid is None or spl_doid == "":
            disease_ancestors_by_spl[spl_name] = {}
        else:
            disease = ontology.search(id="{}".format(spl_doid))[0]
            disease_ancestors_by_spl[spl_name] = {elt.name.replace("DOID_", "DOID:") for elt in disease.ancestors() if elt != owl.Thing}
    return disease_ancestors_by_spl


def getEvidences(record, evidences_by_gene_id, log, annot_field="ANN", assembly="GRCh38"):
    """
    Return list of evidences corresponding to the variant.

    :param record: The annotated variant.
    :type record: anacore.vcf.VCFRecord
    :param evidences_by_gene_id: By gene ENTREZ id the list of evidences.
    :type evidences_by_gene_id: dict
    :param log: Logger of the script.
    :type log: logging.Logger
    :param annot_field: Field used to store annotations in VCF.
    :type annot_field: str
    :param assembly: Assembly used in alignement step. This information is used to check if an imprecise evidence overlap the variant.
    :type assembly: str
    :return: list of evidences corresponding to the variant.
    :rtype: list
    """
    variant_category = getVariantCategory(record)
    processed_associations = set()
    retained_evidences_id = set()
    retained_evidences = list()
    for annot in record.info[annot_field]:
        if annot["HGVSp"] is not None and not annot["HGVSp"].endswith("p.?") and not annot["HGVSp"].replace(")", "").endswith("="):  #################### = exist if variant is locaterd on splice site or promoter
            annot_association_id = "{}:\t{}".format(annot["SYMBOL"], annot["HGVSp"].split(":", 1)[1])
            if annot_association_id not in processed_associations:  # Skip association on other annotation with the same impact
                processed_associations.add(annot_association_id)
                associations = evidences_by_gene_id[annot["Gene"]] if annot["Gene"] in evidences_by_gene_id else []
                log.debug("{}: {} initial associations".format(record.getName(), len(associations)))
                for asso_idx, curr_asso in enumerate(associations):
                    if asso_idx not in retained_evidences_id:
                        if curr_asso["HGVSp_change"] == "" and curr_asso["HGVSp_change_trunc"] == "":  # Imprecise variant: without HGVS
                            if varCatAreCompatible(variant_category, curr_asso["category"]):
                                if curr_asso[assembly + "_category_intervals"] != "":  # Location can be validated
                                    evid_intervals = RegionList()
                                    evid_intervals_str = curr_asso[assembly + "_category_intervals"].split(",")
                                    for curr_interval in evid_intervals_str:
                                        evid_intervals.append(Region.fromStr(curr_interval))
                                    variant_region = Region(
                                        start=int(record.refStart()),  # Extend insertion at the two border nt
                                        end=int(record.refEnd() + 0.5),  # Extend insertion at the two border nt
                                        reference=record.chrom
                                    )  ############################# Manage upstream and downstream
                                    if len(evid_intervals.getOverlapped(variant_region)) > 0:  # Evidence and variant have an overlap
                                        retained_evidences.append(curr_asso)
                                        retained_evidences_id.add(asso_idx)
                                else:
                                    retained_evidences.append(curr_asso)
                                    retained_evidences_id.add(asso_idx)
                        else:  # Precise variant: standard or truncated HGVS
                            record_hgvs_p = HGVSProtChange.fromStr(annot["HGVSp"].rsplit("p.", 1)[1])
                            record_hgvs_p.predicted = False
                            if str(record_hgvs_p) == curr_asso["HGVSp_change"]:  # Standard HGVS with exact match
                                retained_evidences.append(curr_asso)
                                retained_evidences_id.add(asso_idx)
                            elif curr_asso["HGVSp_change"] != "":  ################################# to test
                                # Standard HGVS try to match insertions and duplications/repeats
                                asso_hgvs_p = HGVSProtChange.fromStr(curr_asso["HGVSp_change"])
                                if record_hgvs_p.isInFrameIns() and asso_hgvs_p.isInFrameIns():  # Are insertions
                                    if (record_hgvs_p.evt == "ins" and asso_hgvs_p.evt != "ins"):  # Record is insertion and association is repeat
                                        if HGVSProtChange.insCouldBeIdentical(record_hgvs_p, asso_hgvs_p):
                                            retained_evidences.append(curr_asso)
                                            retained_evidences_id.add(curr_asso["id"])
                                    elif (record_hgvs_p.evt != "ins" and asso_hgvs_p.evt == "ins"):  # Record is repeat and association is insertion
                                        if HGVSProtChange.insCouldBeIdentical(asso_hgvs_p, record_hgvs_p):
                                            retained_evidences.append(curr_asso)
                                            retained_evidences_id.add(asso_idx)
                            else:  # Truncated HGVS (ex: G12)
                                if str(record_hgvs_p).startswith(curr_asso["HGVSp_change_trunc"]):
                                    retained_evidences.append(curr_asso)
                                    retained_evidences_id.add(asso_idx)
    return retained_evidences


def addHeaderFields(writer, disease_id=None, disease_term=None):
    """
    Add evidence fields and disease information in VCF header.

    :param writer: Handle to the VCF output.
    :type writer: anacore.vcf.VCFIO
    :param disease_id: The disease ID in ontology. If this parameter is set, disease ID is apply on all samples.
    :type disease_id: str
    :param disease_term: The disease term or exact synonym in ontology. If this parameter is set, disease term is apply on all samples.
    :type disease_term: str
    """
    writer.format["EVID_PS"] = HeaderFormatAttr("EVID_PS", "The higher clinical evidence level for this precise variant in sample disease.", number=1)
    writer.format["EVID_PA"] = HeaderFormatAttr("EVID_PA", "The higher clinical evidence level for this precise variant in all diseases.", number=1)
    writer.format["EVID_IS"] = HeaderFormatAttr("EVID_IS", "The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in sample disease.", number=1)
    writer.format["EVID_IA"] = HeaderFormatAttr("EVID_IA", "The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in all diseases.", number=1)
    # DOID and DO_term for each sample
    for spl_name in writer.samples:
        if spl_name in writer.sample_info:  # VCF already contains samples information
            spl_attr = writer.sample_info[spl_name]
            # Disease ID
            doid = ""
            if "doid" in spl_attr.keys():
                doid = spl_attr.doid
            else:
                spl_attr.case_by_attr["doid"] = "DOID"
            if disease_id:
                doid = disease_id  # Force disease ID
            spl_attr.doid = doid
            # Disease term
            do_term = ""
            if "do_term" in spl_attr.keys():
                do_term = spl_attr.do_term
            else:
                spl_attr.case_by_attr["do_term"] = "DO_term"
            if disease_term:
                do_term = disease_term  # Force disease term
            spl_attr.do_term = do_term
        else:
            writer.sample_info[spl_name] = HeaderSampleAttr(
                spl_name,
                doid="" if disease_id is None else disease_id,
                doterm="" if disease_term is None else disease_term
            )
            writer.sample_info[spl_name].case_by_attr["doid"] = "DOID"


def writeEvidencesDetails(out_file, details_by_spl):
    """
    Write evidences details on dedicated file.

    :param out_file: Path to the evidences associated with variants. (format: TSV or JSON depends on extension).
    :type out_file: str
    :param details_by_spl: List of evidences by sample name.
    :type details_by_spl: dict
    """
    if not out_file.endswith(".tsv"):
        with open(out_file, "w") as writer_evidences:
            json.dump(details_by_spl, writer_evidences)
    else:
        with HashedSVIO(out_file, "w") as writer_evidences:
            writer_evidences.titles = [
                "sample", "variant_gene", "variant_HGVSp", "variant_genomic_alteration",
                "evidence_subject", "evidence_disease", "evidence_level",
                "evidence_type", "evidence_direction", "evidence_drugs",
                "evidence_clinical_significance", "evidence_citation", "evidence_source",
                "evidence_doid", "evidence_category", "evidence_variant_specificity",
                "evidence_disease_specificity"
            ]
            for spl_name, spl in details_by_spl.items():
                for curr_variant in spl["variants"]:
                    for curr_evidence in curr_variant["evidences"]:
                        curr_evidence["sample"] = spl_name
                        for key in ["gene", "HGVSp", "genomic_alteration"]:
                            curr_evidence["variant_" + key] = curr_variant[key]
                        for key in [elt for elt in writer_evidences.titles if elt.startswith("evidence_") and not elt.endswith("_specificity")]:
                            curr_evidence[key] = curr_evidence[key.replace("evidence_", "")]
                        curr_evidence["evidence_variant_specificity"] = curr_evidence["specificity"]["variant"]
                        curr_evidence["evidence_disease_specificity"] = curr_evidence["specificity"]["disease"]
                        writer_evidences.write(curr_evidence)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Add clinical/drug evidence level to known variants.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-y', '--assembly-version', default="GRCh38", help='Assembly used in alignement step. This information is used to check if an imprecise evidence overlap the variant. [Default: %(default)s]')
    group_disease = parser.add_mutually_exclusive_group()
    group_disease.add_argument('-d', '--disease-id', help='Replace samples disease by this value. It must be a valid ID in disease ontology.')
    group_disease.add_argument('-t', '--disease-term', help='Replace samples disease by this value. It must be a valid TERM in disease ontology.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-e', '--input-evidences', required=True, help='Path to the clinical evidence file from AnaCore-utils/bin/enrichCivic.py (format: TSV).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file. (format: VCF).')
    group_input.add_argument('-n', '--input-disease-ontology', required=True, help='Path to the disease ontology file. (format: OWL).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the annotated file. (format: VCF).')
    group_output.add_argument('-s', '--output-evidences', help='Path to the evidences associated with variants. (format: TSV or JSON depends on extension).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    nb_variants = {"total": 0, "with_asso": 0, "with_precise_asso": 0}
    ontology = get_ontology("file://{}".format(args.input_disease_ontology)).load()
    disease_id, disease_term = getDiseaseElt(ontology, args.disease_id, args.disease_term)
    evidences_by_gene_id = loadEvidences(args.input_evidences)
    details_by_spl = {}
    with AnnotVCFIO(args.output_variants, "w") as writer:
        with AnnotVCFIO(args.input_variants) as reader:
            # Header
            writer.copyHeader(reader)
            addHeaderFields(writer, disease_id, disease_term)
            writer.writeHeader()
            disease_ancestors_by_spl = getDiseaseAncestors(ontology, writer)
            for spl_name, spl in writer.sample_info.items():
                details_by_spl[spl.id] = {
                    "variants": [],
                    "disease": {"doid": spl.doid, "term": spl.doterm}
                }
            del(ontology)  # Ontology is no longer useful
            # Records
            for record in reader:
                log.debug(
                    "start {} ({})".format(
                        record.getName(),
                        sorted(set(annot["SYMBOL"] for annot in record.info[args.annotation_field] if annot["SYMBOL"]))
                    )
                )
                nb_variants["total"] += 1
                # Get corresponding clinical evidences
                retained_evidences = getEvidences(record, evidences_by_gene_id, log, args.annotation_field, args.assembly_version)
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
                        "gene": ";".join(
                            sorted(set(
                                [elt["SYMBOL"] for elt in record.info[args.annotation_field] if elt["SYMBOL"]]
                            ))
                        ),
                        "genomic_alteration": record.getName(),
                        "HGVSp": ";".join(
                            sorted(set(
                                [elt["HGVSp"] for elt in record.info[args.annotation_field] if elt["HGVSp"]]
                            ))
                        ),
                        "evidences": []
                    })
                contains_precise = False
                if len(retained_evidences) != 0:
                    nb_variants["with_asso"] += 1
                    for curr_evidence in sorted(retained_evidences, key=lambda elt: elt["level"]):
                        detailed_annot = {
                            "subject": curr_evidence["subject"],
                            "category": curr_evidence["category"],
                            "disease": curr_evidence['disease'],
                            "doid": curr_evidence["doid"],
                            "level": curr_evidence['level'],
                            "type": curr_evidence['type'],
                            "direction": curr_evidence['direction'],
                            "drugs": curr_evidence['drugs'],
                            "clinical_significance": curr_evidence['clinical_significance'],
                            "citation": curr_evidence['citation'],
                            "source": curr_evidence['source'],
                            "specificity": {"variant": "imprecise", "disease": "all"}
                        }
                        for spl_name in record.samples:
                            disease_ancestors = disease_ancestors_by_spl[spl_name]
                            evidence_for_spl = copy.deepcopy(detailed_annot)
                            # Best evidence
                            variant_specificity = ["i"]
                            if curr_evidence["HGVSp_change"] != "":
                                variant_specificity = ["p", "i"]
                                evidence_for_spl["specificity"]["variant"] = "precise"
                                contains_precise = True
                            disease_specificity = ["a"]
                            if curr_evidence['doid'] in disease_ancestors:
                                disease_specificity = ["a", "s"]
                                evidence_for_spl["specificity"]["disease"] = "specific"
                            for curr_v_spe in variant_specificity:
                                for curr_d_spe in disease_specificity:
                                    best_tag = "{}{}".format(curr_v_spe, curr_d_spe)
                                    prev = best_evidence[spl_name][best_tag]
                                    if prev == "" or prev > curr_evidence['level']:
                                        best_evidence[spl_name][best_tag] = curr_evidence['level']
                            # Add to detail
                            details_by_spl[spl_name]["variants"][-1]["evidences"].append(evidence_for_spl)
                # Add best evidence to record
                for spl_name in record.samples:
                    spl_disease = writer.sample_info[spl_name].doid
                    for key, val in best_evidence[spl_name].items():
                        if key[1] == "s" and spl_disease == "":  # Evidence for disease when sample disease is not provided
                            record.samples[spl_name]["EVID_{}".format(key.upper())] = None
                        else:  # Evidence for all diseases or when sample disease is provided
                            record.samples[spl_name]["EVID_{}".format(key.upper())] = val
                record.format.extend(["EVID_PS", "EVID_PA", "EVID_IS", "EVID_IA"])
                if contains_precise:
                    nb_variants["with_precise_asso"] += 1
                # Write record
                writer.write(record)

    # Write detailed outputs
    if args.output_evidences:
        writeEvidencesDetails(args.output_evidences, details_by_spl)

    # Log
    log.info("{:.2%} of variants have an association ({}/{})".format(
        0 if nb_variants["total"] == 0 else nb_variants["with_asso"] / nb_variants["total"],
        nb_variants["with_asso"],
        nb_variants["total"]
    ))
    log.info("{:.2%} of variants have a precise association ({}/{})".format(
        0 if nb_variants["total"] == 0 else nb_variants["with_precise_asso"] / nb_variants["total"],
        nb_variants["with_precise_asso"],
        nb_variants["total"]
    ))
    log.info("End of job")
