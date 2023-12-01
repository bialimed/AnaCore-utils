#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.3.0'

import os
import sys
import uuid
import pysam
import tempfile
import unittest
from anacore.vcf import VCFRecord, HeaderFormatAttr, HeaderInfoAttr
from anacore.sequenceIO import IdxFastaIO

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from mergeCoOccurVar import getIncludingReadsDNA, getIncludingReadsRNA, \
    getSupportingReads, mergedRecord, setRefPos


########################################################################
#
# FUNCTIONS
#
########################################################################
class LoggerSilencer:
    def debug(self, args):
        pass

    def info(self, args):
        pass


class getIncludingReads(unittest.TestCase):
    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta_path, self.tmp_faidx_path, self.tmp_sam_path, self.tmp_bam_path, self.tmp_bam_path + ".bai"]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_sam_path = os.path.join(tmp_folder, unique_id + ".sam")
        self.tmp_bam_path = os.path.join(tmp_folder, unique_id + ".bam")
        #               11  15   20        30        40        50        60        70        80        90        100
        # pos 1 3 5 7 9 | 13|    |         |         |         |   54    |         |      77 |         |         |
        # ref TCTGGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCATCT
        # aln    GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAA///////////////////////GGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA
        # Write BAM
        with open(self.tmp_sam_path, "w") as FH_sam:
            FH_sam.write("""@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
read1	0	chr1	4	60	50M23N50M	*	0	0	GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAAGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:0	MD:Z:100	AS:i:90	XS:i:0
read2	0	chr1	4	60	123M	*	0	0	GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:0	MD:Z:123	AS:i:113	XS:i:0
read3	0	chr1	4	60	50M23N50M	*	0	0	GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAAGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:0	MD:Z:100	AS:i:90	XS:i:0""")
        with pysam.AlignmentFile(self.tmp_sam_path) as FH_sam:
            with pysam.AlignmentFile(self.tmp_bam_path, "wb", template=FH_sam) as FH_bam:
                for rec in FH_sam:
                    FH_bam.write(rec)
        pysam.index(self.tmp_bam_path)
        # Reference seq
        self.tmp_fasta_path = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_faidx_path = os.path.join(tmp_folder, unique_id + ".fa.fai")
        with open(self.tmp_fasta_path, "w") as FH_seq:
            FH_seq.write(""">chr1
TCTGGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCATCT""")
        with open(self.tmp_faidx_path, "w") as FH_faidx:
            FH_faidx.write("chr1\t123\t6\t200\t201")

    def testRNA(self):
        with pysam.AlignmentFile(self.tmp_bam_path) as FH_aln:
            with IdxFastaIO(self.tmp_fasta_path) as FH_seq:
                # Start of alignments
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 1, 2)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 1, 1)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 2, 2)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 3, 3)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 3, 8)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 4, 4)),
                    3
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 4, 8)),
                    3
                )
                # End of alignments
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 129, 129)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 127, 127)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 125, 128)),
                    0
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 126, 126)),
                    3
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 123, 126)),
                    3
                )
                # Start of skip
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 53, 53)),
                    3
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 54, 54)),
                    1
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 50, 57)),
                    1
                )
                # End of skip
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 76, 76)),
                    1
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 77, 77)),
                    3
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 71, 80)),
                    1
                )
                # Over skip
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 53, 77)),
                    1
                )
                self.assertEqual(
                    len(getIncludingReadsRNA(FH_aln, FH_seq, "chr1", 54, 76)),
                    1
                )


class SetSupportingReads(unittest.TestCase):
    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta_path, self.tmp_faidx_path, self.tmp_sam_path, self.tmp_bam_path, self.tmp_bam_path + ".bai"]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_sam_path = os.path.join(tmp_folder, unique_id + ".sam")
        self.tmp_bam_path = os.path.join(tmp_folder, unique_id + ".bam")
        self.tmp_fasta_path = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_faidx_path = os.path.join(tmp_folder, unique_id + ".fa.fai")
        self.ref_seq = "ggaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgcattggggtg"
        #               | | | | | |  |  |  |  |
        #               1 3 5 7 9 11 14 17 20 23
        with open(self.tmp_fasta_path, "w") as FH_seq:
            FH_seq.write(">chr1\n{}".format(self.ref_seq))
        with open(self.tmp_faidx_path, "w") as FH_faidx:
            FH_faidx.write("chr1\t{}\t6\t200\t201".format(len(self.ref_seq)))
        self.reads_content = """>subtit_AAA/CAC_1_alt
ggaagccctgatcACGCCACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_2_alt
aagccctgatcACGCCACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtit_AAA/CAC_4_mixUp
ggaagccctgatcACGCCAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAA/CAC_5_mixDown
ggaagccctgatcACGCAACTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_1_alt
ggaagccctgatcACGCCCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_2_alt
aagccctgatcACGCCCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtitClose_AA/CC_4_mixUp
ggaagccctgatcACGCCAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtitClose_AA/CC_5_mixDown
ggaagccctgatcACGCACATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_1_alt
ggaagccctgatcACGCCCTTCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_2_alt
aagccctgatcACGCCCTTCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>subtit_AAATCTC/CCTTCGG_4_mixUp
ggaagccctgatcACGCCCTTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>subtit_AAATCTC/CCTTCGG_5_mixDown
ggaagccctgatcACGCAAATCGGGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_1_alt
ggaagccctgatcACGCTGGAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_2_alt
aagccctgatcACGCTGGAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insertion_A/TGGAGG_4_mixUp
ggaagccctgatcACGCTGGAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insertion_A/TGGAGG_5_mixDown
ggaagccctgatcACGCAGGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_1_alt
ggaagccctgatcACGCTGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_2_alt
aagccctgatcACGCTGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>deletion_AAATCTC/T_4_mixUp
ggaagccctgatcACGCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>deletion_AAATCTC/T_5_mixDown
ggaagccctgatcACGCAAATGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_1_alt
ggaagccctgatcACGCTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_2_alt
aagccctgatcACGCTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>delIns_AAAT/TGA_4_mixUp
ggaagccctgatcACGCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delIns_AAAT/TGA_5_mixDown
ggaagccctgatcACGCAAATGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_1_alt
ggaagccctgatcACGCGGGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_2_alt
aagccctgatcACGCGGGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insDel_AAA/GGGA_4_mixUp
ggaagccctgatcACGCGGGAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDel_AAA/GGGA_5_mixDown
ggaagccctgatcACGCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_1_alt
ggaagccctgatcACGCCTGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_2_alt
aagccctgatcACGCCTGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>delInsNoStd_AAATCTC/CTGGG_4_mixUp
ggaagccctgatcACGCCTCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>delInsNoStd_AAATCTC/CTGGG_5_mixDown
ggaagccctgatcACGCAAATGGGCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_AAAT/GTGA_1_alt
ggaagccctgatcACGCGTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_AAAT/GTGA_2_alt
aagccctgatcACGCGTGACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_AAAT/GTGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insDelNoStd_AAAT/GTGA_4_mixUp
ggaagccctgatcACGCGTGAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_AAAT/GTGA_5_mixDown
ggaagccctgatcACGCAACTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_CAAA/CGTGA_1_alt
ggaagccctgatcACGCGTGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_CAAA/CGTGA_2_alt
aagccctgatcACGCGTGATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_CAAA/CGTGA_3_ref
gaagccctgatcACGCAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgc
>insDelNoStd_CAAA/CGTGA_4_mixUp
ggaagccctgatcACGCGTGAAATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca
>insDelNoStd_CAAA/CGTGA_5_mixDown
ggaagccctgatcACGCATCTCGGCATGCCGATTaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcagcaaattcaaccaccagaacattgttcgctgca"""
        self.test_cases = [
            [
                VCFRecord("chr1", 18, "subtit_AAA/CAC", "A", ["C"]),
                VCFRecord("chr1", 20, "subtit_AAA/CAC", "A", ["C"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtit_AAA/CAC_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:17A1A103	AS:i:113	XS:i:0
subtit_AAA/CAC_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:17A105	AS:i:118	XS:i:0
subtit_AAA/CAC_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCAACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:19A103	AS:i:118	XS:i:0
subtit_AAA/CAC_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtit_AAA/CAC_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCACTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:15A1A103	AS:i:111	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "subtitClose_AA/CC", "A", ["C"]),
                VCFRecord("chr1", 19, "subtitClose_AA/CC", "A", ["C"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtitClose_AA/CC_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:17A0A104	AS:i:113	XS:i:0
subtitClose_AA/CC_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:17A105	AS:i:118	XS:i:0
subtitClose_AA/CC_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCACATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:1	MD:Z:18A104	AS:i:118	XS:i:0
subtitClose_AA/CC_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtitClose_AA/CC_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:15A0A104	AS:i:111	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "subtit_AAATCTC/CCTTCGG", "AAA", ["CCT"]),
                VCFRecord("chr1", 23, "subtit_AAATCTC/CCTTCGG", "TC", ["GG"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
subtit_AAATCTC/CCTTCGG_1_alt	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCTTCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17A0A0A2T0C99	AS:i:99	XS:i:0
subtit_AAATCTC/CCTTCGG_4_mixUp	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCCCTTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17A0A0A103	AS:i:108	XS:i:0
subtit_AAATCTC/CCTTCGG_5_mixDown	0	chr1	1	60	123M	*	0	0	GGAAGCCCTGATCACGCAAATCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:22T0C99	AS:i:113	XS:i:0
subtit_AAATCTC/CCTTCGG_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
subtit_AAATCTC/CCTTCGG_2_alt	0	chr1	3	60	121M	*	0	0	AAGCCCTGATCACGCCCTTCGGGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15A0A0A2T0C99	AS:i:99	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "insertion_A/TGGAGG", "-", ["TGG"]),
                VCFRecord("chr1", 19, "insertion_A/TGGAGG", "-", ["GG"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insertion_A/TGGAGG_1_alt	0	chr1	1	60	17M3I1M2I105M	*	0	0	GGAAGCCCTGATCACGCTGGAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:123	AS:i:107	XS:i:0
insertion_A/TGGAGG_4_mixUp	0	chr1	1	60	17M3I106M	*	0	0	GGAAGCCCTGATCACGCTGGAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
insertion_A/TGGAGG_5_mixDown	0	chr1	1	60	18M2I105M	*	0	0	GGAAGCCCTGATCACGCAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:123	AS:i:115	XS:i:0
insertion_A/TGGAGG_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insertion_A/TGGAGG_2_alt	0	chr1	3	60	15M3I1M2I105M	*	0	0	AAGCCCTGATCACGCTGGAGGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:121	AS:i:105	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "deletion_AAATCTC/T", "AAA", ["-"]),
                VCFRecord("chr1", 22, "deletion_AAATCTC/T", "CTC", ["-"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
deletion_AAATCTC/T_1_alt	0	chr1	1	60	17M3D1M3D99M	*	0	0	GGAAGCCCTGATCACGCTGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17^AAA1^CTC99	AS:i:100	XS:i:0
deletion_AAATCTC/T_4_mixUp	0	chr1	1	60	17M3D103M	*	0	0	GGAAGCCCTGATCACGCTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^AAA103	AS:i:111	XS:i:0
deletion_AAATCTC/T_5_mixDown	0	chr1	1	60	21M3D99M	*	0	0	GGAAGCCCTGATCACGCAAATGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:21^CTC99	AS:i:111	XS:i:0
deletion_AAATCTC/T_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
deletion_AAATCTC/T_2_alt	0	chr1	3	60	15M3D1M3D99M	*	0	0	AAGCCCTGATCACGCTGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15^AAA1CTC99	AS:i:99	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "delIns_AAAT/TGA", "AAA", ["-"]),
                VCFRecord("chr1", 22, "delIns_AAAT/TGA", "-", ["GA"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
delIns_AAAT/TGA_1_alt	0	chr1	1	60	17M3D1M2I102M	*	0	0	GGAAGCCCTGATCACGCTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:17^AAA103	AS:i:105	XS:i:0
delIns_AAAT/TGA_4_mixUp	0	chr1	1	60	17M3D103M	*	0	0	GGAAGCCCTGATCACGCTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^AAA103	AS:i:111	XS:i:0
delIns_AAAT/TGA_5_mixDown	0	chr1	1	60	21M2I102M	*	0	0	GGAAGCCCTGATCACGCAAATGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:123	AS:i:115	XS:i:0
delIns_AAAT/TGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
delIns_AAAT/TGA_2_alt	0	chr1	3	60	15M3D1M2I102M	*	0	0	AAGCCCTGATCACGCTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:15^AAA103	AS:i:103	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "insDel_AAA/GGGA", "-", ["GGG"]),
                VCFRecord("chr1", 19, "insDel_AAA/GGGA", "AA", ["-"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insDel_AAA/GGGA_1_alt	0	chr1	1	60	17M3I1M2D103M	*	0	0	GGAAGCCCTGATCACGCGGGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:18^AA103	AS:i:106	XS:i:0
insDel_AAA/GGGA_4_mixUp	0	chr1	1	60	17M3I106M	*	0	0	GGAAGCCCTGATCACGCGGGAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
insDel_AAA/GGGA_5_mixDown	0	chr1	1	60	18M2D103M	*	0	0	GGAAGCCCTGATCACGCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:18^AA103	AS:i:113	XS:i:0
insDel_AAA/GGGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insDel_AAA/GGGA_2_alt	0	chr1	3	60	15M3I1M2D103M	*	0	0	AAGCCCTGATCACGCGGGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:16^AA103	AS:i:104	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "delInsNoStd_AAATCTC/CTGGG", "AAA", ["C"]),
                VCFRecord("chr1", 22, "delInsNoStd_AAATCTC/CTGGG", "-", ["GGG"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
delInsNoStd_AAATCTC/CTGGG_1_alt	0	chr1	1	60	17M2D2M3I102M	*	0	0	GGAAGCCCTGATCACGCCTGGGCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:6	MD:Z:17^AA0C103	AS:i:102	XS:i:0
delInsNoStd_AAATCTC/CTGGG_4_mixUp	0	chr1	1	60	17M2D104M	*	0	0	GGAAGCCCTGATCACGCCTCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^AA0C103	AS:i:108	XS:i:0
delInsNoStd_AAATCTC/CTGGG_5_mixDown	0	chr1	1	60	21M3I102M	*	0	0	GGAAGCCCTGATCACGCAAATGGGCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
delInsNoStd_AAATCTC/CTGGG_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
delInsNoStd_AAATCTC/CTGGG_2_alt	0	chr1	3	60	15M2D2M3I102M	*	0	0	AAGCCCTGATCACGCCTGGGCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:6	MD:Z:15^AA0C103	AS:i:102	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 18, "insDelNoStd_AAAT/GTGA", "A", ["GTG"]),
                VCFRecord("chr1", 20, "insDelNoStd_AAAT/GTGA", "AT", ["-"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insDelNoStd_AAAT/GTGA_1_alt	0	chr1	1	60	17M1D3I1M2D102M	*	0	0	GGAAGCCCTGATCACGCGTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:6	MD:Z:17^A103	AS:i:103	XS:i:0
insDelNoStd_AAAT/GTGA_4_mixUp	0	chr1	1	60	17M1D3I105M	*	0	0	GGAAGCCCTGATCACGCGTGAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:17^A105	AS:i:110	XS:i:0
insDelNoStd_AAAT/GTGA_5_mixDown	0	chr1	1	60	19M2D102M	*	0	0	GGAAGCCCTGATCACGCAACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:19^AT102	AS:i:113	XS:i:0
insDelNoStd_AAAT/GTGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insDelNoStd_AAAT/GTGA_2_alt	0	chr1	3	60	15M1D3I1M2D102M	*	0	0	AAGCCCTGATCACGCGTGACTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:6	MD:Z:15A0A0A0T102	AS:i:102	XS:i:0"""
            ],
            [
                VCFRecord("chr1", 17, "insDelNoStd_CAAA/CGTGA", "C", ["CGTG"]),
                VCFRecord("chr1", 18, "insDelNoStd_CAAA/CGTGA", "AAA", ["A"]),
                """@SQ	SN:chr1	LN:131
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fa reads.fa
insDelNoStd_CAAA/CGTGA_1_alt	0	chr1	1	60	17M3I1M2D103M	*	0	0	GGAAGCCCTGATCACGCGTGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:18^AA103	AS:i:106	XS:i:0
insDelNoStd_CAAA/CGTGA_4_mixUp	0	chr1	1	60	17M3I106M	*	0	0	GGAAGCCCTGATCACGCGTGAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:3	MD:Z:123	AS:i:114	XS:i:0
insDelNoStd_CAAA/CGTGA_5_mixDown	0	chr1	1	60	18M2D103M	*	0	0	GGAAGCCCTGATCACGCATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:2	MD:Z:18^AA103	AS:i:113	XS:i:0
insDelNoStd_CAAA/CGTGA_3_ref	0	chr1	2	60	118M	*	0	0	GAAGCCCTGATCACGCAAATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGC	*	NM:i:0	MD:Z:118	AS:i:118	XS:i:0
insDelNoStd_CAAA/CGTGA_2_alt	0	chr1	3	60	15M3I1M2D103M	*	0	0	AAGCCCTGATCACGCGTGATCTCGGCATGCCGATTAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGATCATCAGCAAATTCAACCACCAGAACATTGTTCGCTGCA	*	NM:i:5	MD:Z:16^AA103	AS:i:104	XS:i:0"""
            ]
        ]

    def testSetSupportingReads(self):
        with IdxFastaIO(self.tmp_fasta_path) as FH_seq:
            for first, second, aln_content in self.test_cases:
                # Write BAM
                with open(self.tmp_sam_path, "w") as FH_sam:
                    FH_sam.write(aln_content)
                with pysam.AlignmentFile(self.tmp_sam_path) as FH_sam:
                    with pysam.AlignmentFile(self.tmp_bam_path, "wb", template=FH_sam) as FH_bam:
                        for rec in FH_sam:
                            FH_bam.write(rec)
                pysam.index(self.tmp_bam_path)
                # Eval
                first.normalizeSingleAllele()
                second.normalizeSingleAllele()
                with pysam.AlignmentFile(self.tmp_bam_path) as FH_aln:
                    setRefPos(first, FH_seq)
                    setRefPos(second, FH_seq)
                    first.isIns = first.isInsertion()
                    second.isIns = second.isInsertion()
                    shared_reads = getIncludingReadsDNA(FH_aln, FH_seq, "chr1", first.upstream_start, second.downstream_end)
                    first.supporting_reads = getSupportingReads(first, FH_seq, "chr1", FH_aln, LoggerSilencer()) & shared_reads
                    second.supporting_reads = getSupportingReads(second, FH_seq, "chr1", FH_aln, LoggerSilencer()) & shared_reads
                    # Check supporting first
                    expected = sorted([
                        "{}_{}".format(first.id, curr_suffix) for curr_suffix in ["1_alt", "2_alt", "4_mixUp"]
                    ])
                    self.assertEqual(
                        sorted(first.supporting_reads),
                        expected
                    )
                    # Check supporting second
                    expected = sorted([
                        "{}_{}".format(second.id, curr_suffix) for curr_suffix in ["1_alt", "2_alt", "5_mixDown"]
                    ])
                    self.assertEqual(
                        sorted(second.supporting_reads),
                        expected
                    )


class FakeVCFIO:
    def __init__(self, info, format):
        self.info = info
        self.format = format


class MergeCoOccurVar(unittest.TestCase):
    def setUp(self):
        # VCF
        self.vcfio = FakeVCFIO(
            {"AF": HeaderInfoAttr("AF", "Alternative alleles frequencies", "Float", "A")},
            {
                "AD": HeaderFormatAttr("AD", "Alternative alleles depths", "Integer", "A"),
                "DP": HeaderFormatAttr("DP", "total depth", "Integer", "1")
            }
        )
        # Ref seq
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_fasta_path = os.path.join(tmp_folder, unique_id + ".fa")
        self.tmp_faidx_path = os.path.join(tmp_folder, unique_id + ".fa.fai")
        self.ref_seq = "ACGCAAATCTCGGCATGCCGATT"
        #               | | | | | |  |  |  |  |
        #               1 3 5 7 9 11 14 17 20 23
        with open(self.tmp_fasta_path, "w") as FH_seq:
            FH_seq.write(">chr1\n{}".format(self.ref_seq))
        with open(self.tmp_faidx_path, "w") as FH_faidx:
            FH_faidx.write("chr1\t{}\t6\t60\t61".format(len(self.ref_seq)))
        # Variants
        self.variant_1 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            "artificial_1",  # id
            None,  # ref
            None,  # alt
            10,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.05]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [10], "DP": 100},
                "splB": {"AD": [40], "DP": 4900},
            }
        )
        self.variant_2 = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            30,  # qual
            ["PASS"],  # filter
            {"AF": [0.06]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )
        self.expected_merge = VCFRecord(
            "chr1",  # chrom
            None,  # pos
            None,  # id
            None,  # ref
            None,  # alt
            20,  # qual
            ["lowQual", "lowDP"],  # filter
            {"AF": [0.06], "MCO_QUAL": [10, 30], "MCO_VAR": ["chr1:5=A/T", "chr1:20=G/C"]},  # info
            ["DP", "AD"],  # format
            {
                "splA": {"AD": [5], "DP": 50},
                "splB": {"AD": [31], "DP": 550},
            }
        )

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_fasta_path, self.tmp_faidx_path]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testMergedRecord(self):
        test_cases = [
            {  # substit and 14 missing ref between
                "prev": {"pos": 5, "ref": "A", "alt": ["T"]},
                "curr": {"pos": 20, "ref": "G", "alt": ["C"]},
                "merge": {"pos": 5, "ref": "AAATCTCGGCATGCCG", "alt": ["TAATCTCGGCATGCCC"]}
            },
            {  # largeSubstit and 1 missing ref between
                "prev": {"pos": 5, "ref": "AAAT", "alt": ["TGCA"]},
                "curr": {"pos": 10, "ref": "TC", "alt": ["GG"]},
                "merge": {"pos": 5, "ref": "AAATCTC", "alt": ["TGCACGG"]}
            },
            {  # largeCloseSubstit
                "prev": {"pos": 5, "ref": "AAAT", "alt": ["TGCA"]},
                "curr": {"pos": 9, "ref": "CT", "alt": ["GG"]},
                "merge": {"pos": 5, "ref": "AAATCT", "alt": ["TGCAGG"]}
            },
            {  # delIns and 1 missing ref between
                "prev": {"pos": 5, "ref": "AAAT", "alt": ["-"]},
                "curr": {"pos": 10, "ref": "-", "alt": ["GGCATCT"]},
                "merge": {"pos": 5, "ref": "AAATC", "alt": ["CGGCATCT"]}
            },

            {  # insDel and 2 missing ref between
                "prev": {"pos": 5, "ref": "-", "alt": ["GTGTG"]},
                "curr": {"pos": 7, "ref": "ATC", "alt": ["-"]},
                "merge": {"pos": 5, "ref": "AAATC", "alt": ["GTGTGAA"]}
            },

            {  # insDel and 1 missing ref between
                "prev": {"pos": 5, "ref": "-", "alt": ["GTGTG"]},
                "curr": {"pos": 6, "ref": "AA", "alt": ["-"]},
                "merge": {"pos": 5, "ref": "AAA", "alt": ["GTGTGA"]}
            },
            {  # closeInsDel
                "prev": {"pos": 5, "ref": "-", "alt": ["GTGTG"]},
                "curr": {"pos": 5, "ref": "AA", "alt": ["-"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["GTGTG"]}
            },
            {  # delIns and 1 missing ref between
                "prev": {"pos": 5, "ref": "AAA", "alt": ["-"]},
                "curr": {"pos": 9, "ref": "-", "alt": ["GGGT"]},
                "merge": {"pos": 5, "ref": "AAAT", "alt": ["TGGGT"]}
            },
            {  # closeDelIns
                "prev": {"pos": 5, "ref": "AAA", "alt": ["-"]},
                "curr": {"pos": 8, "ref": "-", "alt": ["GGGT"]},
                "merge": {"pos": 5, "ref": "AAA", "alt": ["GGGT"]}
            },
            {  # closeSubstit
                "prev": {"pos": 5, "ref": "A", "alt": ["G"]},
                "curr": {"pos": 6, "ref": "A", "alt": ["T"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["GT"]}
            },
            {  # closeInsAfter
                "prev": {"pos": 5, "ref": "A", "alt": ["G"]},
                "curr": {"pos": 6, "ref": "-", "alt": ["T"]},
                "merge": {"pos": 5, "ref": "A", "alt": ["GT"]}
            },
            {  # closeDelAfter
                "prev": {"pos": 5, "ref": "A", "alt": ["G"]},
                "curr": {"pos": 6, "ref": "A", "alt": ["-"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["G"]}
            },
            {  # closeInsBefore
                "prev": {"pos": 6, "ref": "-", "alt": ["T"]},
                "curr": {"pos": 6, "ref": "A", "alt": ["G"]},
                "merge": {"pos": 6, "ref": "A", "alt": ["TG"]}
            },
            {  # closeDelBefore
                "prev": {"pos": 5, "ref": "A", "alt": ["-"]},
                "curr": {"pos": 6, "ref": "A", "alt": ["G"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["G"]}
            },
            {  # delBefore
                "prev": {"pos": 5, "ref": "A", "alt": ["-"]},
                "curr": {"pos": 6, "ref": "A", "alt": ["G"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["G"]}
            },
            {  # coDelIns
                "prev": {"pos": 5, "ref": "AAAT", "alt": ["-"]},
                "curr": {"pos": 9, "ref": "-", "alt": ["AGG"]},
                "merge": {"pos": 5, "ref": "AAAT", "alt": ["AGG"]}
            },
            {  # coInsDel
                "prev": {"pos": 5, "ref": "-", "alt": ["GTGTG"]},
                "curr": {"pos": 5, "ref": "AA", "alt": ["-"]},
                "merge": {"pos": 5, "ref": "AA", "alt": ["GTGTG"]}
            }
        ]
        for curr_test in test_cases:
            self.variant_1.pos = curr_test["prev"]["pos"]
            self.variant_1.ref = curr_test["prev"]["ref"]
            self.variant_1.alt = curr_test["prev"]["alt"]
            self.variant_2.pos = curr_test["curr"]["pos"]
            self.variant_2.ref = curr_test["curr"]["ref"]
            self.variant_2.alt = curr_test["curr"]["alt"]
            self.expected_merge.pos = curr_test["merge"]["pos"]
            self.expected_merge.ref = curr_test["merge"]["ref"]
            self.expected_merge.alt = curr_test["merge"]["alt"]
            self.expected_merge.info = {
                "AF": [0.06],
                "MCO_QUAL": [10, 30],
                "MCO_VAR": [self.variant_1.getName(), self.variant_2.getName()]
            }
            # Eval
            with IdxFastaIO(self.tmp_fasta_path) as FH_ref:
                observed_merge = mergedRecord(
                    self.vcfio,
                    self.variant_1, self.variant_1.getName(),
                    self.variant_2, self.variant_2.getName(),
                    FH_ref)
            self.assertEqual(
                strVariant(observed_merge),
                strVariant(self.expected_merge)
            )


def strVariant(var):
    info = []
    for info_key, info_val in sorted(var.info.items()):
        info.append(
            "{}:{}".format(info_key, info_val)
        )
    samples = []
    for spl_name, spl_val in sorted(var.samples.items()):
        for key, val in sorted(spl_val.items()):
            samples.append(
                "{}:{}:{}".format(spl_name, key, val)
            )
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        var.getName(),
        var.id,
        var.qual,
        sorted(var.filter),
        info,
        sorted(var.format),
        samples
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
