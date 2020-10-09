#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFPrimers(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_noises = os.path.join(tmp_folder, unique_id + ".tsv")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Create VCF
        vcf_content = """##fileformat=VCFv4.1
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	10	allNoise	A	-	.	.	.	DP:AF	500:0.05	500:0.05
1	20	allOver	A	-	.	.	.	DP:AF	500:0.3	500:0.3
1	30	oneOver	A	-	.	.	.	DP:AF	500:0.05	500:0.3
1	40	notKnown	A	-	.	.	.	DP:AF	500:0.05	500:0.05"""
        with open(self.tmp_variants, "w") as FH_variants:
            FH_variants.write(vcf_content)

        # Create TSV
        tsv_content = """#Chromosome	Position	Reference_allele	Alternative_allele	Noise_rate
1	10	A	-	0.05
1	20	A	-	0.05
1	30	A	-	0.05"""
        with open(self.tmp_noises, "w") as FH_noises:
            FH_noises.write(tsv_content)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_noises, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testTag(self):
        cmd = [
            "filterVCFNoise.py",
            "--mode", "tag",
            "--input-variants", self.tmp_variants,
            "--input-noises", self.tmp_noises,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """##fileformat=VCFv4.1
##FILTER=<ID=popConst,Description="The variant correspond to a constitutive detection (sequencing polymerase error rate at this position, workflow artifact, ...).">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	10	allNoise	A	-	.	popConst	.	DP:AF	500:0.05	500:0.05
1	20	allOver	A	-	.	PASS	.	DP:AF	500:0.3	500:0.3
1	30	oneOver	A	-	.	PASS	.	DP:AF	500:0.05	500:0.3
1	40	notKnown	A	-	.	PASS	.	DP:AF	500:0.05	500:0.05""".split("\n")
        observed = list()
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip().split("\n")
        self.assertEqual(
            expected[1:],  # skip file format version
            observed[1:]  # skip file format version
        )

    def testRemove(self):
        cmd = [
            "filterVCFNoise.py",
            "--mode", "remove",
            "--input-variants", self.tmp_variants,
            "--input-noises", self.tmp_noises,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """##fileformat=VCFv4.1
##FILTER=<ID=popConst,Description="The variant correspond to a constitutive detection (sequencing polymerase error rate at this position, workflow artifact, ...).">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	20	allOver	A	-	.	PASS	.	DP:AF	500:0.3	500:0.3
1	30	oneOver	A	-	.	PASS	.	DP:AF	500:0.05	500:0.3
1	40	notKnown	A	-	.	PASS	.	DP:AF	500:0.05	500:0.05""".split("\n")
        observed = list()
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip().split("\n")
        self.assertEqual(
            expected[1:],  # skip file format version
            observed[1:]  # skip file format version
        )

    def testTagChange(self):
        cmd = [
            "filterVCFNoise.py",
            "--mode", "remove",
            "--tag-name", "noise",
            "--tag-description", "The signal is lower than noise.",
            "--input-variants", self.tmp_variants,
            "--input-noises", self.tmp_noises,
            "--output-variants", self.tmp_output
        ]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """##fileformat=VCFv4.1
##FILTER=<ID=noise,Description="The signal is lower than noise.">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA	splB
1	20	allOver	A	-	.	PASS	.	DP:AF	500:0.3	500:0.3
1	30	oneOver	A	-	.	PASS	.	DP:AF	500:0.05	500:0.3
1	40	notKnown	A	-	.	PASS	.	DP:AF	500:0.05	500:0.05""".split("\n")
        observed = list()
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip().split("\n")
        self.assertEqual(
            expected[1:],  # skip file format version
            observed[1:]  # skip file format version
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
