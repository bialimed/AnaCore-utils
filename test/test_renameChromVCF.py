#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
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
class TestRenameChromVCF(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_names = os.path.join(tmp_folder, unique_id + "_names.tsv")
        self.tmp_input = os.path.join(tmp_folder, unique_id + "_in.vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_input, self.tmp_names, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testDefault(self):
        # Create input
        content = """##fileformat=VCFv4.1
##fileDate=20200202
##source=GenerateSVCandidates 1.6.0
##reference=file://Homo_sapiens.GRCh38.94.dna.fa
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
##FILTER=<ID=LowEvidence,Description="RNA fusion calls without both split read and spanning pair support">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
chr1	154944376	1	A	]chr7:72948502]A	.	PASS	SVTYPE=BND;MATEID=10	PR	6,3
chr7	20204330	2	A	]14:20205909]A	.	PASS	SVTYPE=BND;MATEID=11	PR:SR	108,18:14,13
17	20205909	3	C	[chr16:20204330[	.	PASS	SVTYPE=BND;MATEID=12	PR:SR	108,18:14,13
chr20	2020	4	G	T,A	.	PASS	.	AF:DP	0.2,0.1:20"""
        with open(self.tmp_input, "w") as writer:
            writer.write(content)
        # Execute command
        cmd = [
            "renameChromVCF.py",
            "--input-variants", self.tmp_input,
            "--output-variants", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Validate results
        expected = """##fileformat=VCFv4.1
##fileDate=20200202
##source=GenerateSVCandidates 1.6.0
##reference=file://Homo_sapiens.GRCh38.94.dna.fa
##contig=<ID=1,length=248956422>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=2,length=242193529>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=MT,length=16569>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FILTER=<ID=LowEvidence,Description="RNA fusion calls without both split read and spanning pair support">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
1	154944376	1	A	]7:72948502]A	.	PASS	MATEID=10;SVTYPE=BND	PR	6,3
7	20204330	2	A	]14:20205909]A	.	PASS	MATEID=11;SVTYPE=BND	PR:SR	108,18:14,13
17	20205909	3	C	[16:20204330[	.	PASS	MATEID=12;SVTYPE=BND	PR:SR	108,18:14,13
20	2020	4	G	T,A	.	PASS	.	AF:DP	0.2,0.1:20"""
        observed = None
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            expected,
            observed
        )

    def testRenameFile(self):
        # Create rename rules
        content = """1	chr1
10	chr10
11	chr11
12	chr12
13	chr13
14	14
15	chr15
16	chr16
18	chr18
19	chr19
2	chr2
20	chr20
21	chr21
22	chr22
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
M	chrM
X	chrX
Y	chrY"""
        with open(self.tmp_names, "w") as writer:
            writer.write(content)
        # Create input
        content = """##fileformat=VCFv4.1
##fileDate=20200202
##source=GenerateSVCandidates 1.6.0
##reference=file://Homo_sapiens.GRCh38.94.dna.fa
##contig=<ID=1,length=248956422>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=2,length=242193529>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=M,length=16569>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
##FILTER=<ID=LowEvidence,Description="RNA fusion calls without both split read and spanning pair support">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
1	154944376	1	A	]11:72948502]A	.	PASS	SVTYPE=BND;MATEID=10	PR	6,3
11	20204330	2	A	]14:20205909]A	.	PASS	SVTYPE=BND;MATEID=11	PR:SR	108,18:14,13
17	20205909	3	C	[16:20204330[	.	PASS	SVTYPE=BND;MATEID=12	PR:SR	108,18:14,13
20	2020	4	G	T,A	.	PASS	.	AF:DP	0.2,0.1:20"""
        with open(self.tmp_input, "w") as writer:
            writer.write(content)
        # Execute command
        cmd = [
            "renameChromVCF.py",
            "--input-names", self.tmp_names,
            "--input-variants", self.tmp_input,
            "--output-variants", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Validate results
        expected = """##fileformat=VCFv4.1
##fileDate=20200202
##source=GenerateSVCandidates 1.6.0
##reference=file://Homo_sapiens.GRCh38.94.dna.fa
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FILTER=<ID=LowEvidence,Description="RNA fusion calls without both split read and spanning pair support">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
chr1	154944376	1	A	]chr11:72948502]A	.	PASS	MATEID=10;SVTYPE=BND	PR	6,3
chr11	20204330	2	A	]14:20205909]A	.	PASS	MATEID=11;SVTYPE=BND	PR:SR	108,18:14,13
17	20205909	3	C	[chr16:20204330[	.	PASS	MATEID=12;SVTYPE=BND	PR:SR	108,18:14,13
chr20	2020	4	G	T,A	.	PASS	.	AF:DP	0.2,0.1:20"""
        observed = None
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            expected,
            observed
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
