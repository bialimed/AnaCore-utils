#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest
from anacore.vcf import VCFRecord

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(BIN_DIR)

from filterBND import AnnotGetter, inNormal, isHLA, isIG, isInner, isReadthrough, loadNormalDb

os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFilterBND(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_annot = os.path.join(tmp_folder, unique_id + "_annot.gtf")
        self.tmp_normal_db1 = os.path.join(tmp_folder, unique_id + "_normDb1.tsv")
        self.tmp_normal_db2 = os.path.join(tmp_folder, unique_id + "_normDb2.tsv")
        with open(self.tmp_annot, "w") as writer:
            writer.write("""chr1	simulation	exon	100	140	.	+	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "1"; gene_name "GENE_N01";
chr1	simulation	exon	100	140	.	+	.	gene_id "GENE_I04"; transcript_id "TR_04"; exon_number "1"; gene_name "GENE_N04";
chr1	simulation	exon	145	148	.	-	.	gene_id "GENE_I05"; transcript_id "TR_05"; exon_number "1"; gene_name "GENE_N05";
chr1	simulation	exon	150	180	.	+	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "2"; gene_name "GENE_N01";
chr1	simulation	exon	150	180	.	+	.	gene_id "GENE_I04"; transcript_id "TR_04"; exon_number "2"; gene_name "GENE_N04";
chr1	simulation	exon	200	250	.	+	.	gene_id "GENE_I02"; transcript_id "TR_02"; exon_number "1"; gene_name "GENE_N02";
chr1	simulation	exon	200	250	.	+	.	gene_id "GENE_I04"; transcript_id "TR_04"; exon_number "3"; gene_name "GENE_N04";
chr1	simulation	exon	300	350	.	+	.	gene_id "GENE_I03"; transcript_id "TR_03"; exon_number "1"; gene_name "GENE_N03";
chr1	simulation	exon	290	340	.	-	.	gene_id "GENE_I06"; transcript_id "TR_06"; exon_number "1"; gene_name "GENE_N06";""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_annot, self.tmp_normal_db1, self.tmp_normal_db2]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testLoadNormalDb(self):
        with open(self.tmp_normal_db1, "w") as writer:
            writer.write("""GENE_ID01	GENE_ID03
GENE_ID01	GENE_ID04""")
        with open(self.tmp_normal_db2, "w") as writer:
            writer.write("""GENE_ID01	GENE_ID03
GENE_ID02	GENE_ID05""")
        expected = {"GENE_ID01	GENE_ID03", "GENE_ID01	GENE_ID04", "GENE_ID02	GENE_ID05"}
        observed = loadNormalDb([self.tmp_normal_db1, self.tmp_normal_db2])
        self.assertEqual(observed, expected)

    def testInNormal(self):
        normal_fusions = {"GENE_ID01	GENE_ID02", "GENE_ID02	GENE_ID03"}
        up = VCFRecord(
            "chr1", 140, "id_01", "A", ["A[chr1:199["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "ANN": [
                    {"SYMBOL": "GENE_N01", "Gene": "GENE_ID01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "Gene": "GENE_ID04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 199, "id_02", "A", ["]chr1:140]A"],
            info={
                "MATEID": "id_01",
                "ANN": [
                    {"SYMBOL": "GENE_N02", "Gene": "GENE_ID02", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(inNormal(up, down, "ANN", normal_fusions))
        up = VCFRecord(
            "chr1", 140, "id_01", "A", ["A[chr1:299["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "ANN": [
                    {"SYMBOL": "GENE_N01", "Gene": "GENE_ID01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "Gene": "GENE_ID04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 299, "id_02", "A", ["]chr1:140]A"],
            info={
                "MATEID": "id_01",
                "ANN": [
                    {"SYMBOL": "GENE_N03", "Gene": "GENE_ID03", "STRAND": "1"},
                    {"SYMBOL": "GENE_N06", "Gene": "GENE_ID06", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(not inNormal(up, down, "ANN", normal_fusions))
        down = VCFRecord(
            "chr1", 140, "id_01", "A", ["]chr1:299]A"],
            info={
                "MATEID": "id_02",
                "ANN": [
                    {"SYMBOL": "GENE_N01", "Gene": "GENE_ID01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "Gene": "GENE_ID04", "STRAND": "1"}
                ]
            }
        )
        up = VCFRecord(
            "chr1", 299, "id_02", "A", ["A[chr1:140["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_01",
                "ANN": [
                    {"SYMBOL": "GENE_N03", "Gene": "GENE_ID03", "STRAND": "1"},
                    {"SYMBOL": "GENE_N06", "Gene": "GENE_ID06", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(not inNormal(up, down, "ANN", normal_fusions))

    def testInner(self):
        up = VCFRecord(
            "chr1", 140, "id_01", "A", ["A[chr1:299["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 299, "id_02", "A", ["]chr1:140]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N03", "STRAND": "1"},
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(not isInner(up, down, "TESTANN"))  # +/+ not inner (starts on limit)
        up = VCFRecord(
            "chr1", 140, "id_01", "A", ["A[chr1:199["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 199, "id_02", "A", ["]chr1:140]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isInner(up, down, "TESTANN"))  # +/+ inner gene 4 (starts on limit)
        up = VCFRecord(
            "chr1", 298, "id_01", "A", ["]chr1:320]A"],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 320, "id_02", "A", ["A[chr1:298["],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"},
                    {"SYMBOL": "GENE_N03", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isInner(up, down, "TESTANN"))  # -/- inner gene 6
        up = VCFRecord(
            "chr1", 298, "id_01", "A", ["A[chr1:320["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 320, "id_02", "A", ["]chr1:298]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"},
                    {"SYMBOL": "GENE_N03", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(not isInner(up, down, "TESTANN"))  # +/+ inner gene 6 => not valid strand
        up = VCFRecord(
            "chr1", 298, "id_01", "A", ["A[chr1:320["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 320, "id_02", "A", ["A[chr1:298["],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"},
                    {"SYMBOL": "GENE_N03", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isInner(up, down, "TESTANN"))  # +/- inner gene 6

    def testIsIG(self):
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "IGH", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isIG(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "IGH", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": []
            }
        )
        self.assertTrue(isIG(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": []
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": []
            }
        )
        self.assertTrue(not isIG(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "IGF2", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(not isIG(up, down, "TESTANN"))


    def testIsHLA(self):
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "HLA-DRB1", "STRAND": "1"},
                    {"SYMBOL": "HLA-DMB", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isHLA(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "HLA-DRB1", "STRAND": "1"},
                    {"SYMBOL": "HLA-DMB", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": []
            }
        )
        self.assertTrue(isHLA(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": []
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": []
            }
        )
        self.assertTrue(not isHLA(up, down, "TESTANN"))
        up = VCFRecord(
            "chr1", 110, "id_01", "A", ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1", 200, "id_02", "A", ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "HLAN02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(not isHLA(up, down, "TESTANN"))

    def testIsReadthrough(self):
        genes = AnnotGetter(self.tmp_annot)
        up = VCFRecord(
            "chr1",
            110,
            "id_01",
            "A",
            ["A[chr1:200["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_02",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            200,
            "id_02",
            "A",
            ["]chr1:110]A"],
            info={
                "MATEID": "id_01",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isReadthrough(up, down, "TESTANN", genes, 1000))
        up = VCFRecord(
            "chr1",
            110,
            "id_03",
            "A",
            ["A[chr1:300["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_04",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            300,
            "id_04",
            "A",
            ["]chr1:110]A"],
            info={
                "MATEID": "id_03",
                "TESTANN": [
                    {"SYMBOL": "GENE_N03", "STRAND": "1"},
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(not isReadthrough(up, down, "TESTANN", genes, 1000))
        up = VCFRecord(
            "chr1",
            140,
            "id_05",
            "A",
            ["A[chr1:199["],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_06",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            199,
            "id_06",
            "A",
            ["]chr1:140]A"],
            info={
                "MATEID": "id_05",
                "TESTANN": [
                    {"SYMBOL": "GENE_N02", "STRAND": "1"}
                ]
            }
        )
        self.assertTrue(isReadthrough(up, down, "TESTANN", genes, 1000))
        up = VCFRecord(
            "chr1",
            289,
            "id_07",
            "A",
            ["]chr1:148]A"],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_08",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            148,
            "id_08",
            "A",
            ["A[chr1:289["],
            info={
                "MATEID": "id_07",
                "TESTANN": [
                    {"SYMBOL": "GENE_N04", "STRAND": "1"},
                    {"SYMBOL": "GENE_N05", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(isReadthrough(up, down, "TESTANN", genes, 1000))
        up = VCFRecord(
            "chr1",
            180,
            "id_09",
            "A",
            ["]chr1:299]A"],
            info={
                "RNA_FIRST": True,
                "MATEID": "id_10",
                "TESTANN": [
                    {"SYMBOL": "GENE_N01", "STRAND": "1"},
                    {"SYMBOL": "GENE_N04", "STRAND": "1"},
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            299,
            "id_10",
            "A",
            ["A[chr1:180["],
            info={
                "MATEID": "id_09",
                "TESTANN": [
                    {"SYMBOL": "GENE_N03", "STRAND": "1"},
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(not isReadthrough(up, down, "TESTANN", genes, 1000))
        up = VCFRecord(
            "chr1",
            285,
            "id_11",
            "A",
            ["]chr1:148]A"],
            info={
                "RNA_FIRST": True,
                "CIPOS": [0, 4],
                "MATEID": "id_12",
                "TESTANN": [
                    {"SYMBOL": "GENE_N06", "STRAND": "-1"}
                ]
            }
        )
        down = VCFRecord(
            "chr1",
            148,
            "id_12",
            "A",
            ["A[chr1:285["],
            info={
                "MATEID": "id_11",
                "TESTANN": [
                    {"SYMBOL": "GENE_N04", "STRAND": "1"},
                    {"SYMBOL": "GENE_N05", "STRAND": "-1"}
                ]
            }
        )
        self.assertTrue(isReadthrough(up, down, "TESTANN", genes, 1000))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
