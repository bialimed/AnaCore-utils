#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.hgvs import AA_THREE_BY_ONE, HGVSProtChange, ONE_LETTER_AA_LEXIC
from anacore.sv import HashedSVIO
import argparse
import logging
import os
from owlready2 import get_ontology
import re
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAliasesByGeneId(gene_info_path):
    aliases_by_id = {}
    with HashedSVIO(gene_info_path) as reader:
        for record in reader:
            aliases_by_id[record["GeneID"]] = {elt.split(":")[-1] for elt in record["Synonyms"].split("|")}
            aliases_by_id[record["GeneID"]].add(record["Symbol"])
    return aliases_by_id


def cleanedAlias(variant, record, aliases_by_gene_id):
    match = re.search(r"(\w+)-(\w+)", variant)
    if match:  # Fusion
        # variant=BCR-ABL with gene=ABL1
        # variant=MLL-MLLT3 with gene=KMT2A
        first, second = match.groups()
        gene_aliases = aliases_by_gene_id[record["entrez_id"]] if record["entrez_id"] in aliases_by_gene_id else {}
        if first in gene_aliases:
            variant = variant.replace(first + "-", record["gene"] + "-")
        elif second in gene_aliases:
            variant = variant.replace("-" + second, "-" + record["gene"])
    else:  # Other
        gene_aliases = aliases_by_gene_id[record["entrez_id"]] if record["entrez_id"] in aliases_by_gene_id else {}
        words = variant.split(" ")
        for idx, curr_word in enumerate(words):
            if curr_word in gene_aliases and curr_word != record["gene"]:
                words[idx] = record["gene"]
        variant = " ".join(words)
    return variant


def splitMulti(variant):
    variant = variant.strip()
    sub_variants = [variant]
    if " and " in variant:  # , and + doe not exist in CIVIC
        sub_variants = []
        for elt in variant.split(" and "):
            sub_variants.extend(splitMulti(elt))
    elif re.match(r"^" + TRUNC_HGVS_P_CHANGE_REGEXP + "/" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", variant):  # G12/G13
        sub_variants = []
        for elt in variant.split("/"):
            sub_variants.extend(splitMulti(elt))
    elif re.match(r"^\w+-\w+ .+$", variant.lower()) and (HGVSProtChange.isValid(variant.split(" ")[-1]) or re.match(r"^" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", variant)):  # fusion + mut
        sub_variants = []
        for elt in variant.split(" "):
            sub_variants.extend(splitMulti(elt))
    elif re.match(r"^\w+ fusions? .+$", variant.lower()) and (HGVSProtChange.isValid(variant.split(" ")[-1]) or re.match(r"^" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", variant)):  # fusion + mut
        # ALK FUSION I1171
        sub_variants = []
        for elt in variant.rsplit(" ", 1):
            sub_variants.extend(splitMulti(elt))
    elif re.match(r"^[\*\w]\d+[\w\*](/[\w\*])+$", variant):  # Q157P/R
        ref = re.match(r"^([\*\w]\d+)[\w\*](/[\w\*])+$", variant).groups()[0]
        alternatives = variant[len(ref):]
        sub_variants = [ref + curr_alt for curr_alt in alternatives.split("/")]
    return sub_variants


def getCategory(variant, record):
    category = getSimpleCategory(variant, record)
    if category == "mutation" and "c." in variant:
        hgvs_c = [elt for elt in variant.split(" ") if "c." in elt]
        if len(hgvs_c) != 1:
            raise Exception('The HGVSc cannot be extract from "{}" in variant {}.'.format(variant, record["variant_id"]))
        hgvs_c = hgvs_c[0]
        if ("ins" in hgvs_c or "dup" in hgvs_c) and "del" not in hgvs_c:
            category = "mutation>insertion"
        elif "del" in hgvs_c and "ins" not in hgvs_c and "dup" not in hgvs_c:
            category = "mutation>deletion"
        else:
            # (c.272_273delinsAA)
            # (c.532C>T)
            # (c.463+8C>T)
            match = re.search(r"^\(?c\.\d+([\+-]\d+)?([ATGC]+)>([ATGC]+)\)?$", hgvs_c)
            if match:
                offset, ref, alt = match.groups()
                if len(ref) == len(alt):
                    if len(ref) == 1:
                        category = "mutation>snv"
                    else:
                        category = "mutation>mnv"
                else:
                    if len(ref) < len(alt):
                        category = "mutation>insertion"
                    else:
                        category = "mutation>deletion"
    return category


def getSimpleCategory(variant, record):
    category = None
    lc_variant = variant.lower()
    if "loh" == lc_variant or "loss of heterozygosity" in lc_variant:
        # LOH
        # COPY-NEUTRAL LOSS OF HETEROZYGOSITY
        category = "mutation>deletion"
    if "expression" in lc_variant:
        # UNDEREXPRESSION
        # EXPRESSION
        # OVEREXPRESSION
        # p16 EXPRESSION
        # ISOFORM EXPRESSION
        # NUCLEAR EXPRESSION
        if "overexpression" in lc_variant:
            category = "expression>overexpression"
        elif "underexpression" in lc_variant:
            category = "expression>underexpression"
        else:
            category = "expression"
    elif "methylation" in lc_variant:
        # PROMOTER METHYLATION
        # PROMOTER DEMETHYLATION
        # PROMOTER HYPERMETHYLATION
        if "hypermethylation" in lc_variant:
            category = "methylation>hypermethylation"
        elif "demethylation" in lc_variant:
            category = "methylation>demethylation"
        else:
            category = "methylation"
    elif "phosphorylation" in lc_variant:
        # T172 PHOSPHORYLATION
        # PHOSPHORYLATION
        category = "phosphorylation"
    elif "of function" in lc_variant.replace("-", " "):
        # Gain-of-Function
        # LOSS-OF-FUNCTION
        # Loss-of-function
        if lc_variant.startswith("gain") or lc_variant.startswith("loss"):
            category = "mutation"
        else:
            raise Exception('The feature "{}" in variant {} cannot be categorized.'.format(variant, record["variant_id"]))
    elif "copy number" in lc_variant or "copy_number" in lc_variant or "loss" in lc_variant or "gain" in lc_variant or "amplification" in lc_variant:
        # copy number gain
        # LOSS
        # copy number loss
        # AMPLIFICATION
        # COPY NUMBER VARIATION
        if "gain" in lc_variant:
            category = "copy_number>gain"
        elif "amplification" in lc_variant:
            category = "copy_number>gain"
        elif "loss" in lc_variant:
            category = "copy_number>loss"
        else:
            category = "copy_number"
    elif "fusion" in lc_variant:
        # ALK FUSIONS
        # FUSION
        category = "fusion"
    elif (lc_variant.startswith(record["gene"].lower() + "-") or lc_variant.endswith("-" + record["gene"].lower())) and re.match(r"^\w+-\w+$", lc_variant):
        # BCR-ABL
        category = "fusion"
    elif lc_variant == "rearrangement":
        # REARRANGEMENT
        category = "fusion"
    elif re.match(r"^rs\d+$", lc_variant.split(" ")[0]):
        # rs681673
        # RS67376798 HOMOZYGOSITY
        category = "mutation"
        if record["reference_bases"] != "" or record["variant_bases"] != "":
            if len(record["reference_bases"]) == len(record["variant_bases"]):
                if len(record["reference_bases"]) == 1:
                    category = "mutation>substitution"
                else:
                    category = "mutation>substitution"
            else:
                if len(record["reference_bases"]) < len(record["variant_bases"]):
                    category = "mutation>insertion"
                else:
                    category = "mutation>deletion"
    elif lc_variant == "itd" or lc_variant == "internal duplication" or "tandem repeat" in lc_variant or "insertion" in lc_variant.split(" "):
        # ITD
        # INTERNAL DUPLICATION
        # Alu insertion
        # 5' TANDEM REPEAT
        category = "mutation>insertion"
    elif re.match("^" + TRUNC_HGVS_P_CHANGE_REGEXP.lower() + "ins$", lc_variant):
        # Which amino acids are inserted ?
        # P780INS
        # C77ins
        # R108ins
        # R177ins
        category = "mutation>insertion"
    elif "complete deletion" in lc_variant:
        # Null (Complete deletion)
        category = "copy_number>loss"
    elif "deletion" in lc_variant:
        # Null (Partial deletion of Exons 2 & 3)
        # DELETION
        # DELETION POLYMORPHISM
        # EXON 11-19 DELETION
        # EXON 19 DELETION
        # EXON 12-22 DELETION
        # 3' EXON DELETION
        # DELETION
        # IKZF1 deletion and mutation
        category = "mutation>deletion"
    elif "del" in lc_variant.split(" "):
        # Ex19 del L858R
        # DEL 485-490
        category = "mutation>deletion"
    elif "frameshift" in lc_variant.replace("frame shift", "frameshift"):
        # EXON 9 FRAMESHIFT
        # N-TERMINAL FRAME SHIFT
        # FRAMESHIFT TRUNCATION
        # FRAMESHIFT MUTATION
        category = "mutation"  # The frameshift most often come from indel but this is not a rule (change in STS and in splicing).
    elif "wildtype" in lc_variant.replace(" ", ""):
        # WILD TYPE
        # WILDTYPE
        category = "wild type"
    elif HGVSProtChange.isValid(variant) or re.match(r"^" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", variant):
        # D842V
        # R882
        # R233*
        # I774FS
        # *757L
        # E55=
        # V2288FS*1
        # Q56_V60del
        # M1?
        # D842_I843delinsVM
        # L755_T759del
        category = "mutation"
        if record["reference_bases"] != "" or record["variant_bases"] != "":
            if len(record["reference_bases"]) == len(record["variant_bases"]):
                if len(record["reference_bases"]) == 1:
                    category = "mutation>substitution"
                else:
                    category = "mutation>substitution"
            else:
                if len(record["reference_bases"]) < len(record["variant_bases"]):
                    category = "mutation>insertion"
                else:
                    category = "mutation>deletion"
        elif re.match(r"^" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", variant):
            if ("ins" in lc_variant or "dup" in lc_variant) and "del" not in lc_variant:
                category = "mutation>insertion"
            elif "[" in lc_variant:  # Repeated sequences
                category = "mutation>insertion"
            elif "del" in lc_variant and "ins" not in lc_variant and "dup" not in lc_variant:
                category = "mutation>deletion"
        else:
            hgvs_p = HGVSProtChange.fromStr(variant)
            if hgvs_p.evt == "del":
                category = "mutation>deletion"
            elif hgvs_p.evt in {"ins", "del"} or "[" in hgvs_p.new_elts:
                category = "mutation>insertion"
    elif "mutation" in lc_variant or "mut" in lc_variant.split(" ") or lc_variant.startswith("splicing alteration (c."):
        # Splicing alteration (c.340+1G>A)
        # MUTATION
        # FGFR2 mutations
        # PROMOTER MUTATION
        # 3' UTR MUTATION
        # DNA BINDING DOMAIN MUTATION
        # Non-P-Loop Mutation
        # P-Loop Mutation
        # TKD MUTATION
        # SH2 DOMAIN MUTATION
        # KINASE DOMAIN MUTATION
        # B2 DOMAIN MUTATION
        # CONSERVED DOMAIN MUT
        # exon 3 mutations
        # EXON 10 and EXON 21 MUTATIONS
        # EXON 20 INSERTION
        # EXON 1-2 MUTATION
        # RARE EGRF MUT
        # DELETERIOUS MUTATION
        # INACTIVATING MUTATION           # => loss of function
        # TRUNCATING MUTATION
        category = "mutation"
        if record["reference_bases"] != "" or record["variant_bases"] != "":
            if len(record["reference_bases"]) == len(record["variant_bases"]):
                if len(record["reference_bases"]) == 1:
                    category = "mutation>substitution"
                else:
                    category = "mutation>substitution"
            else:
                if len(record["reference_bases"]) < len(record["variant_bases"]):
                    category = "mutation>insertion"
                else:
                    category = "mutation>deletion"
    elif "polymorphism" in lc_variant:
        # 3' UTR Polymorphism
        if record["reference_bases"] != "" or record["variant_bases"] != "":
            if len(record["reference_bases"]) == len(record["variant_bases"]):
                if len(record["reference_bases"]) == 1:
                    category = "mutation>substitution"
                else:
                    category = "mutation>substitution"
            else:
                if len(record["reference_bases"]) < len(record["variant_bases"]):
                    category = "mutation>insertion"
                else:
                    category = "mutation>deletion"
    return category


def getOntologyElementByID(ontology, disease_id):
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


def cleanedRawVariant(record):
    feature_change_by_variant_id = {
        "1991": "EXON 10 MUTATIONS and EXON 21 MUTATIONS",  # Instead of "EXON 10 and EXON 21 MUTATIONS"
        "785": "Thr367Metfs",  # Instead of "1100DELC"
        "736": "V769_D770insASV",  # Instead of "V769_770insASV"
        "503": "EML4-ALK",  # Instead of "EML4-ALK E6;A20"
        "501": "EML4-ALK",  # Instead of "EML4-ALK E2;A20"
        "500": "EML4-ALK",  # Instead of "EML4-ALK E20;A20"
        "352": "EML4-ALK and C1156Y and L1198F",
        "312": "EXON 2-7 DELETION (EGFRvIII)",
        "173": "EML4-ALK T1151dup",  # Instead of "EML4-ALK T1151INST"
        "164": "EWSR1-FLI1",  # Instead of "EWSR1-FLI1 Type 1"
    }
    cleaned_variant = record["variant"]
    if record["variant_id"] in feature_change_by_variant_id:
        cleaned_variant = feature_change_by_variant_id[record["variant_id"]]
    return cleaned_variant


HGVS_P_CHANGE_REGEXP = r"(([\w\*]\d+=)|([\w\*]\d+\w*\*?))"
TRUNC_HGVS_P_CHANGE_REGEXP = r"[" + "".join(ONE_LETTER_AA_LEXIC).replace("*", r"\*") + r"]\d+"


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Standardize clinical evidence from CIViC bank.')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-bank', required=True, help='Path to the clinical evidence file from CIViC (format: TSV).')
    group_input.add_argument('-d', '--input-disease-ontology', required=True, help='Path to the disease ontology file. (format: OWL).')
    group_input.add_argument('-g', '--input-genes', help='Path to the genes information from NCBI-Entrez (format: TSV).')
    group_output = parser.add_argument_group('Outputs')
    group_input.add_argument('-o', '--output-bank', required=True, help='Path to the outputted file (format: TSV).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    renamed_diseases = set()
    ontology = get_ontology("file://{}".format(args.input_disease_ontology)).load()
    aliases_by_gene_id = {}
    if args.input_genes:
        aliases_by_gene_id = getAliasesByGeneId(args.input_genes)
    with HashedSVIO(args.output_bank, "w") as writer:
        writer.titles = ["gene", "entrez_id", "subject", "category", "disease",
                         "doid", "level", "type", "drugs", "direction",
                         "clinical_significance", "citation", "HGVSp_change",
                         "HGVSp_change_trunc", "source", "transcript", "chromosome",
                         "start", "end", "xref_level"]
        with HashedSVIO(args.input_bank) as reader:
            for record in reader:
                if record["evidence_status"] == "accepted":
                    cleaned_variant = cleanedRawVariant(record).replace("  ", " ")
                    splitted_feature = splitMulti(cleaned_variant)
                    # Disease
                    record["disease"] = record["disease"]
                    if record["disease"] == "Pediatric Low-grade Glioma (PLGG)":
                        record["disease"] = "Pediatric Low-grade Glioma"
                    ontology_elt = getOntologyElementByTerm(ontology, record["disease"], False)
                    if ontology_elt is not None:
                        elt_ids = set([str(ontology_elt.id[0])] + [str(elt) for elt in ontology_elt.hasAlternativeId])
                        if record["doid"] != "" and "DOID:" + record["doid"] not in elt_ids:
                            raise Exception('The disease "{}" match to the ID {} in ontology and is associated with DOID:{}.'.format(record["disease"], ontology_elt.id, record["doid"]))
                        record["doid"] = ontology_elt.id[0][5:]
                        record["disease"] = str(ontology_elt.label[0])
                    elif record["doid"] != "":
                        ontology_elt = getOntologyElementByID(ontology, "DOID:{}".format(record["doid"]))
                        if ontology_elt is None:
                            raise Exception("The variant {} does not have valid DO_term and DO_id.".format(record["variant_id"]))
                        disease_term = str(ontology_elt.label[0])
                        if record["disease"] not in renamed_diseases:
                            log.warning('Rename disease "{}" to "{}".'.format(record["disease"], disease_term))
                            renamed_diseases.add(record["disease"])
                        record["disease"] = disease_term
                    else:
                        log.warning('The variant {} with evidence level {} does not have DO_id and DO_term is not found: {}.'.format(record["variant_id"], record["evidence_level"], record["disease"]))
                    # Other fields
                    for sub_variant in splitted_feature:
                        # Clean
                        cleaned_sub_variant = cleanedAlias(sub_variant, record, aliases_by_gene_id)
                        if cleaned_sub_variant.startswith(record["gene"] + " "):
                            cleaned_sub_variant = cleaned_sub_variant[len(record["gene"]) + 1:]
                        if re.match(r"^.+ \(?c\..+\)?$", sub_variant):
                            probable_hgvs = cleaned_sub_variant.split(" ")[0]
                            if HGVSProtChange.isValid(probable_hgvs):
                                cleaned_sub_variant = str(HGVSProtChange.fromStr(probable_hgvs))
                            elif re.match(TRUNC_HGVS_P_CHANGE_REGEXP, cleaned_sub_variant):
                                cleaned_sub_variant = probable_hgvs
                        record["subject"] = cleaned_sub_variant
                        # HGVS
                        record["HGVSp_change"] = ""
                        record["HGVSp_change_trunc"] = ""
                        if HGVSProtChange.isValid(cleaned_sub_variant):
                            record["HGVSp_change"] = HGVSProtChange.fromStr(cleaned_sub_variant)
                        elif re.match("^" + TRUNC_HGVS_P_CHANGE_REGEXP + "$", cleaned_sub_variant):
                            record["HGVSp_change_trunc"] = AA_THREE_BY_ONE[cleaned_sub_variant[0]] + cleaned_sub_variant[1:]
                        # Category
                        category = getCategory(cleaned_sub_variant, record)
                        record["category"] = "" if category is None else category
                        # Position
                        record["transcript"] = record['representative_transcript'] + ("" if record['representative_transcript2'] == "" else " and " + record['representative_transcript2'])
                        record["chromosome"] = ""
                        record["start"] = ""
                        record["end"] = ""
                        # Level
                        record["xref_level"] = record["evidence_level"]
                        record["level"] = record["evidence_level"]
                        if record["evidence_level"] == "C":
                            record["level"] = "C" if record["evidence_type"] == "predictive" else "D"
                        elif record["evidence_level"] == "E":
                            record["level"] = "D"
                        # Additionnal features
                        record["source"] = "CIViC"
                        record["type"] = record["evidence_type"]
                        record["direction"] = record["evidence_direction"]
                        record["citation"] = "{}:{}".format(record['source_type'].replace("PubMed", "PMID"), record['citation_id'])
                        if record['drug_interaction_type'] != "":
                            record["drugs"] = "{} of {}".format(record['drug_interaction_type'], record['drugs'].replace(",", ", "))
                        # Write
                        writer.write(record)
    log.info("End of job")


########################################################################
#
# TO MANAGE
#
########################################################################
# ALTERATION  # P53 ALTERATION is a bucket type variant used in studies which is constructed by combining all cases of p53 mutation together with p53 overexpression
# intron 10 rearrangement  # no fusion
# SNP309
# c.128-?_250+?
# BIALLELIC INACTIVATION
# c.1641+1dup
# CYTOPLASMIC MISLOCALIZATION
# SPLICE VARIANT 7
# NUCLEAR TRANSLOCATION
# ALTERNATIVE TRANSCRIPT (ATI)
# SERUM LEVELS
# Non-V600
# Ex19 del L858R  # pb to retrieve HGVS
# EXON 14 SKIPPING MUTATION # Alternative splicing

# E746_T751>I  # is delins
# P772_V774insPHV   # is delins not ins

# 106insR (c.316insGCC)  # Missing starting amino acids
# 77insL(c.230insTCT)  # Missing starting amino acids and space before "(c."
