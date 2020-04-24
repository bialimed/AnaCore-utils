#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import uuid
import tempfile
import unittest
from anacore.vcf import VCFRecord
from anacore.gtf import loadModel
from anacore.region import splittedByRef

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(BIN_DIR)

from annotBND import annotGeneShard, exonsPos, getMostSupported, shardIsBeforeBND

os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestAnnotBND(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_annot = os.path.join(tmp_folder, unique_id + "_annot.gtf")
        with open(self.tmp_annot, "w") as writer:
            """
            Model:
                              10     40       80 90    140  150   170 180 190      220 230    250
                              |      |        |  |      |   |       |  |  |         |  |      |
            gene 2 tr 3                          >>>>>>>>                              >>>>>>>>
            gene 2 tr 4                          >>>>>>>>>>>>          >>>>            >>>>>>>>
                           --------------------------------------------------------------------
            gene 1 tr 1       <<<<<<<<        <<<<<<<<<<<<<<<       <<<<<<<<<<<<<<<<<
            gene 1 tr 2       <<<<<<<<                              <<<<<<<<<<<<<<<<<
            """
            writer.write("""1	simulation	exon	10	40	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "3"; gene_name "GENE_N01";
1	simulation	exon	10	40	.	-	.	gene_id "GENE_I01"; transcript_id "TR_02"; exon_number "2"; gene_name "GENE_N01";
1	simulation	exon	80	150	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "2"; gene_name "GENE_N01";
1	simulation	exon	90	140	.	+	.	gene_id "GENE_I02"; transcript_id "TR_03"; exon_number "1"; gene_name "GENE_N02";
1	simulation	exon	90	150	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "1"; gene_name "GENE_N02";
1	simulation	exon	170	230	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "1"; gene_name "GENE_N01";
1	simulation	exon	170	230	.	-	.	gene_id "GENE_I01"; transcript_id "TR_02"; exon_number "1"; gene_name "GENE_N01";
1	simulation	exon	180	190	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "2"; gene_name "GENE_N02";
1	simulation	exon	220	250	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "2"; gene_name "GENE_N02";
1	simulation	exon	220	250	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "3"; gene_name "GENE_N02";""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_annot]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_shardIsBeforeBND(self):
        record = VCFRecord(
            "1", 70, "id_01", "A", ["A[2:100["],
            info={"RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertTrue(shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_02", "A", ["A]2:100]"],
            info={"RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertTrue(shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_03", "A", ["[2:100[A"],
            info={"RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertTrue(not shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_04", "A", ["]2:100]A"],
            info={"RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertTrue(not shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_05", "A", ["A[2:100["],
            info={"MATEID": "id_02"}
        )
        self.assertTrue(shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_06", "A", ["A]2:100]"],
            info={"MATEID": "id_02"}
        )
        self.assertTrue(shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_07", "A", ["[2:100[A"],
            info={"MATEID": "id_02"}
        )
        self.assertTrue(not shardIsBeforeBND(record))
        record = VCFRecord(
            "1", 70, "id_08", "A", ["]2:100]A"],
            info={"MATEID": "id_02"}
        )
        self.assertTrue(not shardIsBeforeBND(record))

    def test_annotGeneShard(self):
        # First +/+
        record = VCFRecord(
            "1", 70, "id_01", "A", ["A[2:100["],
            info={"RNA_FIRST": True, "MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'up'}, {"STRAND": '-', "GENE_SHARD": 'down'}]
        )
        # First +/-
        record = VCFRecord(
            "1", 70, "id_02", "A", ["A]2:100]"],
            info={"RNA_FIRST": True, "MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'up'}, {"STRAND": '-', "GENE_SHARD": 'down'}]
        )
        # First -/-
        record = VCFRecord(
            "1", 70, "id_03", "A", ["]2:100]A"],
            info={"RNA_FIRST": True, "MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'down'}, {"STRAND": '-', "GENE_SHARD": 'up'}]
        )
        # First -/+
        record = VCFRecord(
            "1", 70, "id_04", "A", ["[2:100[A"],
            info={"RNA_FIRST": True, "MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'down'}, {"STRAND": '-', "GENE_SHARD": 'up'}]
        )
        # Second +/+
        record = VCFRecord(
            "1", 70, "id_05", "A", ["]2:100]A"],
            info={"MATEID": "id_05", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'down'}, {"STRAND": '-', "GENE_SHARD": 'up'}]
        )
        # Second +/-
        record = VCFRecord(
            "1", 70, "id_06", "A", ["A]2:100]"],
            info={"MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'up'}, {"STRAND": '-', "GENE_SHARD": 'down'}]
        )
        # Second -/-
        record = VCFRecord(
            "1", 70, "id_07", "A", ["A[2:100["],
            info={"MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'up'}, {"STRAND": '-', "GENE_SHARD": 'down'}]
        )
        # Second -/+
        record = VCFRecord(
            "1", 70, "id_08", "A", ["[2:100[A"],
            info={"MATEID": "id_02", "ANN": [{"STRAND": '+'}, {"STRAND": '-'}]}
        )
        annotGeneShard(record, "ANN")
        self.assertEqual(
            record.info["ANN"],
            [{"STRAND": '+', "GENE_SHARD": 'down'}, {"STRAND": '-', "GENE_SHARD": 'up'}]
        )

    def test_exonsPos(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        # strand +, not movable
        record = VCFRecord(
            "1", 70, "id_01", "A", ["A[2:100["],
            info={"RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {})
        # strand +, movable to left and right, corresponding to exons break at right
        record = VCFRecord(
            "1", 70, "id_03", "A", ["A[2:100["],
            info={"CIPOS": [-30, 40], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {90: 2})
        # strand -, movable to left and right, corresponding to exons break at left and right
        record = VCFRecord(
            "1", 70, "id_04", "A", ["]2:100]A"],
            info={"CIPOS": [-30, 40], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {40: 2, 80: 1})
        # strand -, movable to left and right, just befrore exons breaks
        record = VCFRecord(
            "1", 70, "id_05", "A", ["]2:100]A"],
            info={"CIPOS": [-29, 9], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {})
        # strand -, movable to left and right, on exons breaks
        record = VCFRecord(
            "1", 70, "id_06", "A", ["]2:100]A"],
            info={"CIPOS": [-30, 10], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {40: 2, 80: 1})
        # strand +, movable to left and right, just befrore exons breaks
        record = VCFRecord(
            "1", 70, "id_07", "A", ["A[2:100["],
            info={"CIPOS": [-29, 19], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {})
        # strand +, movable to left and right, on exon break
        record = VCFRecord(
            "1", 70, "id_08", "A", ["A[2:100["],
            info={"CIPOS": [-30, 20], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(exonsPos(record, genes_by_chr), {90: 2})

    # def test_getGeneAnnot(self):
    #     genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))

    def test_getMostSupported(self):
        self.assertEqual(getMostSupported({}), None)
        self.assertEqual(getMostSupported({40: 2, 80: 1}), 40)
        # Select most upstream on ambiguous
        self.assertEqual(getMostSupported({80: 2, 40: 2}), 40)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
