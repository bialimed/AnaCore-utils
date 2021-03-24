#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from filterDetailedEvidences import filterJSON, filterTSV, getKeptVariants


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFilterDetailedEvidences(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testFilterJSON(self):
        in_content = '{"splA": {"variants": [{"gene": "TP53", "genomic_alteration": "17:7676051=G/C", "HGVSp": "NP_000537.3:p.(Ser106Arg)", "evidences": [{"subject": "MUTATION", "category": "mutation", "disease": "chronic lymphocytic leukemia", "doid": "DOID:1040", "level": "A", "type": "Prognostic", "direction": "Supports", "drugs": "", "clinical_significance": "Poor Outcome", "citation": "PMID:26837699", "source": "CIViC", "specificity": {"variant": "imprecise", "disease": "all"}}, {"subject": "DELETERIOUS MUTATION", "category": "mutation", "disease": "head and neck squamous cell carcinoma", "doid": "DOID:5520", "level": "B", "type": "Prognostic", "direction": "Supports", "drugs": "", "clinical_significance": "Poor Outcome", "citation": "PMID:22090360", "source": "CIViC", "specificity": {"variant": "imprecise", "disease": "all"}}]}, {"gene": "TERT", "genomic_alteration": "5:1295113=G/A", "HGVSp": "NP_937983.2:p.(=)", "evidences": []}], "disease": {"doid": "", "term": ""}}}'
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        # Filter one of two
        kept_variants_by_spl = {"splA": {"5:1295113=G/A"}, "splB": {"17:7676051=G/C"}}
        expected = '{"splA": {"variants": [{"gene": "TERT", "genomic_alteration": "5:1295113=G/A", "HGVSp": "NP_937983.2:p.(=)", "evidences": []}], "disease": {"doid": "", "term": ""}}}'
        filterJSON(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = reader.readline()
        self.assertEqual(observed, expected)
        # Filter all variants
        kept_variants_by_spl = {"splA": []}
        expected = '{"splA": {"variants": [], "disease": {"doid": "", "term": ""}}}'
        filterJSON(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = reader.readline()
        self.assertEqual(observed, expected)
        # Filter all
        kept_variants_by_spl = {}
        expected = '{}'
        filterJSON(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = reader.readline()
        self.assertEqual(observed, expected)

    def testFilterTSV(self):
        in_content = """sample\tvariant_gene\tvariant_HGVSp\tvariant_genomic_alteration\tevidence_subject\tevidence_disease\tevidence_level\tevidence_type\tevidence_direction\tevidence_drugs\tevidence_clinical_significance\tevidence_citation\tevidence_source\tevidence_doid\tevidence_category\tevidence_variant_specificity\tevidence_disease_specificity
splA\tTP53\tNP_000537.3:p.(Ser106Arg)\t17:7676051=G/C\tMUTATION\tchronic lymphocytic leukemia\tA\tPrognostic\tSupports\t\tPoor Outcome\tPMID:26837699\tCIViC\tDOID:1040\tmutation\timprecise\tall
splA\tTP53\tNP_000537.3:p.(Ser106Arg)\t17:7676051=G/C\tDELETERIOUS MUTATION\thead and neck squamous cell carcinoma\tB\tPrognostic\tSupports\t\tPoor Outcome\tPMID:22090360\tCIViC\tDOID:5520\tmutation\timprecise\tall"""
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        # No filter and one without details
        kept_variants_by_spl = {"splA": {"5:1295113=G/A", "17:7676051=G/C"}, "splB": {"12:45544384=G/T"}}
        expected = """sample\tvariant_gene\tvariant_HGVSp\tvariant_genomic_alteration\tevidence_subject\tevidence_disease\tevidence_level\tevidence_type\tevidence_direction\tevidence_drugs\tevidence_clinical_significance\tevidence_citation\tevidence_source\tevidence_doid\tevidence_category\tevidence_variant_specificity\tevidence_disease_specificity
splA\tTP53\tNP_000537.3:p.(Ser106Arg)\t17:7676051=G/C\tMUTATION\tchronic lymphocytic leukemia\tA\tPrognostic\tSupports\t\tPoor Outcome\tPMID:26837699\tCIViC\tDOID:1040\tmutation\timprecise\tall
splA\tTP53\tNP_000537.3:p.(Ser106Arg)\t17:7676051=G/C\tDELETERIOUS MUTATION\thead and neck squamous cell carcinoma\tB\tPrognostic\tSupports\t\tPoor Outcome\tPMID:22090360\tCIViC\tDOID:5520\tmutation\timprecise\tall"""
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        filterTSV(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)
        # Filter all variants
        kept_variants_by_spl = {"splA": []}
        expected = "sample\tvariant_gene\tvariant_HGVSp\tvariant_genomic_alteration\tevidence_subject\tevidence_disease\tevidence_level\tevidence_type\tevidence_direction\tevidence_drugs\tevidence_clinical_significance\tevidence_citation\tevidence_source\tevidence_doid\tevidence_category\tevidence_variant_specificity\tevidence_disease_specificity"
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        filterTSV(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)
        # Filter all
        kept_variants_by_spl = {}
        expected = "sample\tvariant_gene\tvariant_HGVSp\tvariant_genomic_alteration\tevidence_subject\tevidence_disease\tevidence_level\tevidence_type\tevidence_direction\tevidence_drugs\tevidence_clinical_significance\tevidence_citation\tevidence_source\tevidence_doid\tevidence_category\tevidence_variant_specificity\tevidence_disease_specificity"
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        filterTSV(self.tmp_in, self.tmp_out, kept_variants_by_spl)
        observed = ""
        with open(self.tmp_out, "r") as reader:
            observed = "".join(reader.readlines()).strip()
        self.assertEqual(observed, expected)

    def testGetKeptVariants(self):
        in_content = """##fileformat=VCFv4.3
##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|HGVS_OFFSET|HGVSg|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|CADD_PHRED|CADD_RAW|MetaLR_rankscore|VEST4_rankscore|FILTER">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allele Depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=EVID_IA,Number=1,Type=String,Description="The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in any diseases.">
##FORMAT=<ID=EVID_IS,Number=1,Type=String,Description="The higher clinical evidence level for this precise and imprecise variants in same area and with same category (example: deletion in exon 8) in sample disease.">
##FORMAT=<ID=EVID_PA,Number=1,Type=String,Description="The higher clinical evidence level for this precise variant in any diseases.">
##FORMAT=<ID=EVID_PS,Number=1,Type=String,Description="The higher clinical evidence level for this precise variant in sample disease.">
##SAMPLE=<ID=splA,DOID="",Doterm="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
17	7676051	.	G	C	35936.7	PASS	ANN=C|missense_variant|MODERATE|TP53|7157|Transcript|NM_000546.5|protein_coding|4/11||NM_000546.5%3Ac.318C>G|NP_000537.3%3Ap.(Ser106Arg)|520|318|106|S/R|agC/agG|rs1555526581&CM013441&COSV52748282||-1||EntrezGene|||||||G|G|||tolerated(0.41)|benign(0.149)||NC_000017.11%3Ag.7676051G>C||||||||||||||||||uncertain_significance|0&0&1&1&1|1&1&1&1&1|||||||18.41|1.889636|0.99058|0.34659|PASS	DP:AD:EVID_PS:EVID_PA:EVID_IS:EVID_IA	3008:1451:.::.:A
5	1295113	.	G	A	3522.03	CSQ	ANN=A|upstream_gene_variant|MODIFIER|TERT|7015|Transcript|NM_198253.3|protein_coding|||NM_198253.2%3Ac.-124C>T|NP_937983.2%3Ap.(%3D)||||||rs1242535815|45|-1||EntrezGene|||||||G|G||||||NC_000005.10%3Ag.1295113G>A|||||||||||||||||||||||||||8.231|0.675795|||ANN.CSQ	DP:AD:EVID_PS:EVID_PA:EVID_IS:EVID_IA	1439:317:.::.:"""
        # Sample 2 variants
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content)
        expected = {"splA": {"17:7676051=G/C", "5:1295113=G/A"}}
        self.assertEqual(getKeptVariants(self.tmp_in), expected)
        # Sample 1 variant
        in_content_1_var = "\n".join(in_content.split("\n")[:-1])
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content_1_var)
        expected = {"splA": {"17:7676051=G/C"}}
        self.assertEqual(getKeptVariants(self.tmp_in), expected)
        # Sample empty
        in_content_0_var = "\n".join(in_content_1_var.split("\n")[:-1])
        with open(self.tmp_in, "w") as writer:
            writer.write(in_content_0_var)
        expected = {"splA": set()}
        self.assertEqual(getKeptVariants(self.tmp_in), expected)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
