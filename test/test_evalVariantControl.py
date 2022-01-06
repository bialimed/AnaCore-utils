#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import subprocess
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class EvalVariantControl(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_ctrl = os.path.join(tmp_folder, unique_id + "_ctrl.vcf")
        self.tmp_test = os.path.join(tmp_folder, unique_id + "_test.vcf")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.tsv")
        # Write targets
        with open(self.tmp_ctrl, "w") as writer:
            writer.write("""##fileformat=VCFv4.1
##fileDate=20181129
##reference=hg38
##source=Tru-Q 3 (5% Tier) Reference Standard(catalog_ID:HD730)
##INFO=<ID=GI,Number=1,Type=String,Description="Gene ID">
##INFO=<ID=FC,Number=1,Type=String,Description="Functional consequence">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Minor Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Control
7	55174771	rs727504233	AGGAATTAAGAGAAGC	A	.	PASS	GI=EGFR;FC=746_750del	AF	0.042
7	140753336	rs113488022	A	T	.	PASS	GI=BRAF;FC=V600E	AF	0.08
7	140753337	rs121913378	C	T	.	PASS	GI=BRAF;FC=V600M	AF	0.04
12	25225628	rs121913527	C	T	.	PASS	GI=KRAS;FC=A146T	AF	0.05""")
        # Write alignments
        with open(self.tmp_test, "w") as writer:
            writer.write("""##fileformat=VCFv4.1
##INFO=<ID=STATUS,Number=1,Type=String,Description="Status in comparison">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Minor Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Control
7	5455487	.	A	T	.	PASS	STATUS=FP	AF	0.02
#7	55174771	rs727504233	AGGAATTAAGAGAAGC	A	.	PASS	STATUS=FN	AF	0.042
7	140753336	rs113488022	A	T	.	PASS	STATUS=TP	AF	0.06
7	140753337	rs121913378	C	T	.	PASS	STATUS=TP	AF	0.08
12	25225628	rs121913527	C	T	.	PASS	STATUS=TP	AF	0.04""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_ctrl, self.tmp_test, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_onlyExpected(self):
        # Execute command
        cmd = [
            "evalVariantControl.py",
            "--only-expected",
            "--expected-file", self.tmp_ctrl,
            "--detected-file", self.tmp_test,
            "--output-file", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Comparison
        expected = """[Summary]
#Nb_checked	Errors_sum	Errors_ratio_sum	Error_out_of_threshold_(0.2)	TP	FN	FP
4	0.11200	2.45000	4/4	3	1	0

[Details]
#Chr:Pos	Ref/Alt	Expected	Detected	StatusError	Error_ratio	Out_of_threshold
12:25225628	C/T	0.05000	0.04000	TP	0.01000	0.20000	True
7:140753336	A/T	0.08000	0.06000	TP	0.02000	0.25000	True
7:140753337	C/T	0.04000	0.08000	TP	-0.04000	1.00000	True
7:55174772	GGAATTAAGAGAAGC/-	0.04200	0.00000	FN	0.04200	1.00000	True"""
        observed = None
        with open(self.tmp_out) as reader:
            observed = reader.read().strip()
        self.assertEqual(expected, observed)

    def test_all(self):
        # Execute command
        cmd = [
            "evalVariantControl.py",
            "--expected-file", self.tmp_ctrl,
            "--detected-file", self.tmp_test,
            "--output-file", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Comparison
        expected = """[Summary]
#Nb_checked	Errors_sum	Errors_ratio_sum	Error_out_of_threshold_(0.2)	TP	FN	FP
5	0.13200	3.45000	5/5	3	1	1

[Details]
#Chr:Pos	Ref/Alt	Expected	Detected	StatusError	Error_ratio	Out_of_threshold
12:25225628	C/T	0.05000	0.04000	TP	0.01000	0.20000	True
7:140753336	A/T	0.08000	0.06000	TP	0.02000	0.25000	True
7:140753337	C/T	0.04000	0.08000	TP	-0.04000	1.00000	True
7:5455487	A/T	0.00000	0.02000	FP	-0.02000	1.00000	True
7:55174772	GGAATTAAGAGAAGC/-	0.04200	0.00000	FN	0.04200	1.00000	True"""
        observed = None
        with open(self.tmp_out) as reader:
            observed = reader.read().strip()
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
