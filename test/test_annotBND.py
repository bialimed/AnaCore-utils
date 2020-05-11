#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
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

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from annotBND import annotGeneShard, annotModelRetIntron, exonsPos, getDistBeforeCDSForward, getDistBeforeCDSReverse, getGeneAnnot, getMostSupported, selectedPos, shardIsBeforeBND
# todo: annot


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
                              10     40       80 90    140  150   170 180 189      220 230    250
                              |      |        |  |      |   |       |  |  |         |  |      |
            gene 2 tr 3                          >>>>>>>>                           >>>>>>>>>>>
                   pr 3                                ..                           .......
            gene 2 tr 4                          >>>>>>>>>>>>          >>>>         >>>>>>>>>>>
                   pr 4                                ......          ....         .......
                           --------------------------------------------------------------------
            gene 1 tr 1       <<<<<<<<        <<<<<<<<<<<<<<<       <<<<<<<<<<<<<<<<<<<<
                   pr 1       ........        ......
            gene 1 tr 2       <<<<<<<<                              <<<<<<<<<<<<<<<<<<<<
                   pr 2       ........
            """
            writer.write("""1	simulation	exon	10	40	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "3"; gene_name "GENE_N01";
1	simulation	CDS	10	40	.	-	.	gene_id "GENE_I01"; protein_id "PROT_01"; transcript_id "TR_01"; exon_number "3"; gene_name "GENE_N01";
1	simulation	exon	10	40	.	-	.	gene_id "GENE_I01"; transcript_id "TR_02"; exon_number "2"; gene_name "GENE_N01";
1	simulation	CDS	10	40	.	-	.	gene_id "GENE_I01"; protein_id "PROT_02"; transcript_id "TR_02"; exon_number "2"; gene_name "GENE_N01";
1	simulation	exon	80	150	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "2"; gene_name "GENE_N01";
1	simulation	CDS	80	100	.	-	.	gene_id "GENE_I01"; protein_id "PROT_01"; transcript_id "TR_01"; exon_number "2"; gene_name "GENE_N01";
1	simulation	exon	90	140	.	+	.	gene_id "GENE_I02"; transcript_id "TR_03"; exon_number "1"; gene_name "GENE_N02";
1	simulation	CDS	135	140	.	+	.	gene_id "GENE_I02"; protein_id "PROT_03"; transcript_id "TR_03"; exon_number "1"; gene_name "GENE_N02";
1	simulation	exon	90	150	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "1"; gene_name "GENE_N02";
1	simulation	CDS	135	150	.	+	.	gene_id "GENE_I02"; protein_id "PROT_04"; transcript_id "TR_04"; exon_number "1"; gene_name "GENE_N02";
1	simulation	exon	170	230	.	-	.	gene_id "GENE_I01"; transcript_id "TR_01"; exon_number "1"; gene_name "GENE_N01";
1	simulation	exon	170	230	.	-	.	gene_id "GENE_I01"; transcript_id "TR_02"; exon_number "1"; gene_name "GENE_N01";
1	simulation	exon	180	189	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "2"; gene_name "GENE_N02";
1	simulation	CDS	180	189	.	+	.	gene_id "GENE_I02"; protein_id "PROT_04"; transcript_id "TR_04"; exon_number "2"; gene_name "GENE_N02";
1	simulation	exon	220	250	.	+	.	gene_id "GENE_I02"; transcript_id "TR_03"; exon_number "2"; gene_name "GENE_N02";
1	simulation	CDS	220	240	.	+	.	gene_id "GENE_I02"; protein_id "PROT_03"; transcript_id "TR_03"; exon_number "2"; gene_name "GENE_N02";
1	simulation	exon	220	250	.	+	.	gene_id "GENE_I02"; transcript_id "TR_04"; exon_number "3"; gene_name "GENE_N02";
1	simulation	CDS	220	240	.	+	.	gene_id "GENE_I02"; protein_id "PROT_04"; transcript_id "TR_04"; exon_number "3"; gene_name "GENE_N02";
2	simulation	exon	50	100	.	-	.	gene_id "GENE_I03"; transcript_id "TR_05"; exon_number "2"; gene_name "GENE_N03";
2	simulation	exon	150	200	.	-	.	gene_id "GENE_I03"; transcript_id "TR_05"; exon_number "1"; gene_name "GENE_N03";""")

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


    def test_getDistBeforeCDSReverse(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        protein = genes_by_chr["1"][0].children[1].proteins[0]
        self.assertEqual(getDistBeforeCDSReverse(230, protein), 111)
        self.assertEqual(getDistBeforeCDSReverse(150, protein), 50)
        self.assertEqual(getDistBeforeCDSReverse(100, protein), 0)

    def test_getDistBeforeCDSForward(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        protein = genes_by_chr["1"][1].children[0].proteins[0]
        self.assertEqual(getDistBeforeCDSForward(90, protein), 45)
        self.assertEqual(getDistBeforeCDSForward(135, protein), 0)

    def test_selectedPos(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        # strand +/+ without exon boundary overlap
        first = VCFRecord(
            "1", 295, "id_01", "A", ["[1:395[A"],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 395, "id_02", "A", ["A]1:295]"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (295, 395))
        # strand -/-
        first = VCFRecord(
            "2", 95, "id_01", "A", ["[1:145[A"],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 145, "id_02", "A", ["A]2:95]"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (100, 150))
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        # strand +/+
        first = VCFRecord(
            "1", 135, "id_01", "A", ["A[1:215["],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 215, "id_02", "A", ["[1:135[A"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (140, 220))
        # strand +/-
        first = VCFRecord(
            "1", 135, "id_01", "A", ["A]1:45]"],
            info={"CIPOS": [0, 10], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 35, "id_02", "A", ["A]1:145]"],
            info={"CIPOS": [0, 10], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (140, 40))
        # strand +/- no common => arbitrary selection of first
        first = VCFRecord(
            "1", 135, "id_01", "A", ["A]1:41]"],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 35, "id_02", "A", ["A]1:141]"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (140, 36))
        # strand +/- only first has exon boundary in cipos
        first = VCFRecord(
            "1", 135, "id_01", "A", ["A]1:401]"],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 395, "id_02", "A", ["A]1:141]"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (140, 396))
        # strand +/+ only first has exon boundary in cipos
        first = VCFRecord(
            "1", 135, "id_01", "A", ["A[1:395["],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 395, "id_02", "A", ["]1:135]A"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (140, 400))
        # strand +/- only second has exon boundary in cipos
        first = VCFRecord(
            "1", 5, "id_01", "A", ["A[1:41["],
            info={"CIPOS": [0, 6], "RNA_FIRST": True, "MATEID": "id_02"}
        )
        second = VCFRecord(
            "1", 35, "id_02", "A", ["A]1:11]"],
            info={"CIPOS": [0, 6], "MATEID": "id_01"}
        )
        observed = selectedPos(
            first, exonsPos(first, genes_by_chr),
            second, exonsPos(second, genes_by_chr)
        )
        self.assertEqual(observed, (6, 40))

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

    def test_annotModelRetIntron(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        # Fragments +/-
        # tr3->tr1&tr2: intron&CDS->spliceAcceptor&CDS undetermined
        # tr4->tr1&tr2: spliceDonor&CDS->spliceAcceptor&CDS frameshift
        record = VCFRecord("1", 189, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 189, "RNA_FIRST": True})
        mate = VCFRecord("1", 40, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 40})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0'),
                ('TR_02', 'TR_01:0&TR_02:0'),
                ('TR_03', 'TR_01:.&TR_02:.'),
                ('TR_04', 'TR_01:0&TR_02:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:.&TR_04:0'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:.&TR_04:0')
            ]
        )
        # Fragments +/-
        # tr3->tr1&tr2: spliceDonor&CDS->spliceAcceptor&CDS in_frame
        # tr4->tr1&tr2: CDS->spliceAcceptor&CDS in_frame
        record = VCFRecord("1", 140, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 140, "RNA_FIRST": True})
        mate = VCFRecord("1", 40, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 40})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0'),
                ('TR_02', 'TR_01:0&TR_02:0'),
                ('TR_03', 'TR_01:1&TR_02:1'),
                ('TR_04', 'TR_01:1&TR_02:1')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:1&TR_04:1'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:1&TR_04:1')
            ]
        )
        # Fragments +/-
        # tr3->tr1&tr2: intron&CDS->CDS undetermined, intron&CDS->intron&UTR in_frame
        # tr4->tr1&tr2: CDS->CDS in_frame, CDS->intron&UTR undetermined
        record = VCFRecord("1", 188, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 188, "RNA_FIRST": True})
        mate = VCFRecord("1", 84, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 84})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0'),
                ('TR_02', 'TR_01:0&TR_02:0'),
                ('TR_03', 'TR_01:.&TR_02:1'),
                ('TR_04', 'TR_01:1&TR_02:.')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:.&TR_04:1'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:1&TR_04:.')
            ]
        )
        # Fragments +/-
        # tr3->tr1&tr2: intron&CDS->intron&CDS in_frame, intron&CDS->intron&UTR in_frame
        # tr4->tr1&tr2: intron&CDS->intron&CDS frameshift, intron&CDS->intron&UTR frameshift
        record = VCFRecord("1", 195, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 195, "RNA_FIRST": True})
        mate = VCFRecord("1", 60, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 60})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0'),
                ('TR_02', 'TR_01:0&TR_02:0'),
                ('TR_03', 'TR_01:1&TR_02:1'),
                ('TR_04', 'TR_01:0&TR_02:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:1&TR_04:0'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:1&TR_04:0')
            ]
        )
        # Fragments +/-
        # tr3->tr1&tr2: spliceDonor&CDS->spliceAcceptor&UTR in_frame, spliceDonor&CDS->spliceAcceptor&UTR frameshift
        # tr4->tr1&tr2: CDS->spliceAcceptor&UTR in_frame, CDS->spliceAcceptor&UTR frameshift
        record = VCFRecord("1", 140, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 140, "RNA_FIRST": True})
        mate = VCFRecord("1", 230, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 230})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:0&TR_04:0'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:0&TR_04:0'),
                ('TR_03', 'TR_01:1&TR_02:0&TR_03:0&TR_04:0'),
                ('TR_04', 'TR_01:1&TR_02:0&TR_03:0&TR_04:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_01:0&TR_02:0&TR_03:1&TR_04:1'),
                ('TR_02', 'TR_01:0&TR_02:0&TR_03:0&TR_04:0'),
                ('TR_03', 'TR_01:0&TR_02:0&TR_03:0&TR_04:0'),
                ('TR_04', 'TR_01:0&TR_02:0&TR_03:0&TR_04:0')
            ]
        )
        # Fragments +/-
        # tr3->tr1&tr2: 3'UTR->CDS undetermined, 3'UTR->intron undetermined
        # tr4->tr1&tr2: 3'UTR->CDS undetermined, 3'UTR->intron undetermined
        record = VCFRecord("1", 245, "id_01", "A", ["A]1:40]"], info={"ANNOT_POS": 245, "RNA_FIRST": True})
        mate = VCFRecord("1", 83, "id_02", "A", ["A[1:189["], info={"ANNOT_POS": 83})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_03', 'TR_01:.&TR_02:.'),
                ('TR_04', 'TR_01:.&TR_02:.')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_03:.&TR_04:.'),
                ('TR_02', 'TR_03:.&TR_04:.')
            ]
        )
        # Fragments -/-
        # tr5->tr1&tr2: intron&untranslated->intron&5'UTR in_frame
        record = VCFRecord("2", 120, "id_01", "A", ["]1:160]A"], info={"ANNOT_POS": 120, "RNA_FIRST": True})
        mate = VCFRecord("1", 160, "id_02", "A", ["A]2:120]"], info={"ANNOT_POS": 160})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_05', 'TR_01:1&TR_02:1&TR_03:0&TR_04:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_05:1'),
                ('TR_02', 'TR_05:1'),
                ('TR_03', 'TR_05:0'),
                ('TR_04', 'TR_05:0')
            ]
        )
        # Fragments -/-
        # tr5->tr1&tr2: exon&spliceDonor&untranslated->exon&spliceAcceptor&5'UTR in_frame,
        #               exon&spliceDonor&untranslated->intron&5'UTR undetermined
        record = VCFRecord("2", 150, "id_01", "A", ["]1:150]A"], info={"ANNOT_POS": 150, "RNA_FIRST": True})
        mate = VCFRecord("1", 150, "id_02", "A", ["A]2:150]"], info={"ANNOT_POS": 150})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_05', 'TR_01:1&TR_02:.&TR_03:0&TR_04:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_05:1'),
                ('TR_02', 'TR_05:.'),
                ('TR_03', 'TR_05:0'),
                ('TR_04', 'TR_05:0')
            ]
        )
        # Fragments -/-
        # tr5->tr1&tr2: exon&spliceDonor&untranslated->exon&5'UTR undetermined,
        #               exon&spliceDonor&untranslated->exon&5'UTR undetermined
        record = VCFRecord("2", 150, "id_01", "A", ["]1:190]A"], info={"ANNOT_POS": 150, "RNA_FIRST": True})
        mate = VCFRecord("1", 190, "id_02", "A", ["A]2:150]"], info={"ANNOT_POS": 190})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_05', 'TR_01:.&TR_02:.&TR_03:0&TR_04:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_01', 'TR_05:.'),
                ('TR_02', 'TR_05:.'),
                ('TR_03', 'TR_05:0'),
                ('TR_04', 'TR_05:0')
            ]
        )
        # Fragments -/+
        # tr5->tr1&tr2: exon&spliceDonor&untranslated->exon&3'UTR no_frame,
        #               exon&spliceDonor&untranslated->exon&3'UTR no_frame
        record = VCFRecord("2", 150, "id_01", "A", ["]1:241]A"], info={"ANNOT_POS": 150, "RNA_FIRST": True})
        mate = VCFRecord("1", 241, "id_02", "A", ["]2:150]A"], info={"ANNOT_POS": 241})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_05', 'TR_03:0&TR_04:0')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_03', 'TR_05:0'),
                ('TR_04', 'TR_05:0')
            ]
        )
        # Fragments -/+
        # tr5->tr1&tr2: exon&spliceDonor&untranslated->exon&CDS undetermined,
        #               exon&spliceDonor&untranslated->exon&CDS undetermined
        record = VCFRecord("2", 150, "id_01", "A", ["]1:235]A"], info={"ANNOT_POS": 150, "RNA_FIRST": True})
        mate = VCFRecord("1", 235, "id_02", "A", ["]2:150]A"], info={"ANNOT_POS": 235})
        record.info["ANN"] = sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"])
        mate.info["ANN"] = sorted(getGeneAnnot(mate, genes_by_chr), key=lambda elt: elt["Feature"])
        annotModelRetIntron(record, mate, "ANN")
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in record.info["ANN"]],  # observed
            [  # expected
                ('TR_05', 'TR_03:.&TR_04:.')
            ]
        )
        self.assertEqual(
            [(elt["Feature"], elt["IN_FRAME"]) for elt in mate.info["ANN"]],  # observed
            [  # expected
                ('TR_03', 'TR_05:.'),
                ('TR_04', 'TR_05:.')
            ]
        )


    def test_getGeneAnnot(self):
        genes_by_chr = splittedByRef(loadModel(self.tmp_annot, "genes"))
        # Strand +, intergenic, before breakend
        record = VCFRecord(
            "1", 5, "id_01", "A", ["A[2:100["],
            info={"ANNOT_POS": 5, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            getGeneAnnot(record, genes_by_chr),
            []
        )
        # Strand -, intergenic, after breakend
        record = VCFRecord(
            "1", 5, "id_02", "A", ["]2:100]A"],
            info={"ANNOT_POS": 5, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            getGeneAnnot(record, genes_by_chr),
            []
        )
        # Strand +, intergenic, after breakend
        record = VCFRecord(
            "1", 5, "id_03", "A", ["]2:100]A"],
            info={"ANNOT_POS": 5, "MATEID": "id_02"}
        )
        self.assertEqual(
            getGeneAnnot(record, genes_by_chr),
            []
        )
        # Strand -, intergenic, before breakend
        record = VCFRecord(
            "1", 5, "id_04", "A", ["A[2:100["],
            info={"ANNOT_POS": 5, "MATEID": "id_02"}
        )
        self.assertEqual(
            getGeneAnnot(record, genes_by_chr),
            []
        )
        # Strand +, intron, before breakend
        record = VCFRecord(
            "1", 50, "id_05", "A", ["A[2:100["],
            info={"ANNOT_POS": 50, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )
        # Strand -, intron, after breakend
        record = VCFRecord(
            "1", 50, "id_06", "A", ["]2:100]A"],
            info={"ANNOT_POS": 50, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 21,
                    "Protein_position": 7,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron&utr",
                    "RNA_ELT_POS": "1/1&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 0
                }
            ]
        )
        # Strand +, intron, after breakend
        record = VCFRecord(
            "1", 50, "id_07", "A", ["]2:100]A"],
            info={"ANNOT_POS": 50, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )
        # Strand -, intron, before breakend
        record = VCFRecord(
            "1", 50, "id_08", "A", ["A[2:100["],
            info={"ANNOT_POS": 50, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 22,
                    "Protein_position": 8,
                    "Codon_position": 1
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron&utr",
                    "RNA_ELT_POS": "1/1&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 0
                }
            ]
        )
        # Strand -, exon spliceDonnor, after breakend
        record = VCFRecord(
            "1", 80, "id_09", "A", ["]2:100]A"],
            info={"ANNOT_POS": 80, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&spliceDonor",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 21,
                    "Protein_position": 7,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron&utr",
                    "RNA_ELT_POS": "1/1&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 0
                }
            ]
        )
        # Strand -, exon spliceAcceptor, after breakend
        record = VCFRecord(
            "1", 40, "id_10", "A", ["]2:100]A"],
            info={"ANNOT_POS": 40, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&spliceAcceptor",
                    "RNA_ELT_POS": "3/3",
                    "CDS_position": 22,
                    "Protein_position": 8,
                    "Codon_position": 1
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&spliceAcceptor",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 1,
                    "Protein_position": 1,
                    "Codon_position": 1
                }
            ]
        )
        # Strand +, exon spliceDonnor, before breakend
        record = VCFRecord(
            "1", 189, "id_11", "A", ["A[2:100["],
            info={"ANNOT_POS": 189, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 70
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 20
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": 6,
                    "Protein_position": 2,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&spliceDonor",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 26,
                    "Protein_position": 9,
                    "Codon_position": 2
                }
            ]
        )
        # Strand +, exon spliceAcceptor, after breakend
        record = VCFRecord(
            "1", 220, "id_12", "A", ["]2:100]A"],
            info={"ANNOT_POS": 220, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 101
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 51
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&spliceAcceptor",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 7,
                    "Protein_position": 3,
                    "Codon_position": 1
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&spliceAcceptor",
                    "RNA_ELT_POS": "3/3",
                    "CDS_position": 27,
                    "Protein_position": 9,
                    "Codon_position": 3
                }
            ]
        )
        # Strand +, exon transcriptStart on strand -, after breakend
        record = VCFRecord(
            "1", 230, "id_13", "A", ["]2:100]A"],
            info={"ANNOT_POS": 230, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&transcriptStart&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 111
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&transcriptStart&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 61
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 17,
                    "Protein_position": 6,
                    "Codon_position": 2
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "3/3",
                    "CDS_position": 37,
                    "Protein_position": 13,
                    "Codon_position": 1
                }
            ]
        )
        # Strand +, exon transcriptStart, before breakend
        record = VCFRecord(
            "1", 90, "id_14", "A", ["A]2:100]"],
            info={"ANNOT_POS": 90, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 11,
                    "Protein_position": 4,
                    "Codon_position": 2
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron&utr",
                    "RNA_ELT_POS": "1/1&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 0
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&transcriptStart&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 45
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&transcriptStart&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 45
                }
            ]
        )
        # Strand -, exon, after breakend
        record = VCFRecord(
            "1", 39, "id_15", "A", ["]2:100]A"],
            info={"ANNOT_POS": 39, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "3/3",
                    "CDS_position": 23,
                    "Protein_position": 8,
                    "Codon_position": 2
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/2",
                    "CDS_position": 2,
                    "Protein_position": 1,
                    "Codon_position": 2
                }
            ]
        )
        # Strand -, exon, after breakend
        record = VCFRecord(
            "1", 92, "id_16", "A", ["]2:100]A"],
            info={"ANNOT_POS": 92, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 9,
                    "Protein_position": 3,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron&utr",
                    "RNA_ELT_POS": "1/1&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 0
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 43
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 43
                }
            ]
        )
        # Strand +, exon, before breakend
        record = VCFRecord(
            "1", 92, "id_17", "A", ["A[2:100["],
            info={"ANNOT_POS": 92, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 9,
                    "Protein_position": 3,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 43
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 43
                }
            ]
        )
        # Strand +, exon, before breakend
        record = VCFRecord(
            "1", 138, "id_18", "A", ["A[2:100["],
            info={"ANNOT_POS": 138, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 38
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "1/2",
                    "CDS_position": 4,
                    "Protein_position": 2,
                    "Codon_position": 1
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "1/3",
                    "CDS_position": 4,
                    "Protein_position": 2,
                    "Codon_position": 1
                }
            ]
        )
        # Strand +, exon
        record = VCFRecord(
            "1", 143, "id_19", "A", ["A[2:100["],
            info={"ANNOT_POS": 143, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 43
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": 6,
                    "Protein_position": 2,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "1/3",
                    "CDS_position": 9,
                    "Protein_position": 3,
                    "Codon_position": 3
                }
            ]
        )
        # Strand +, exon, before breakend
        record = VCFRecord(
            "1", 183, "id_20", "A", ["A[2:100["],
            info={"ANNOT_POS": 183, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 64
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 14
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": 6,
                    "Protein_position": 2,
                    "Codon_position": 3
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 20,
                    "Protein_position": 7,
                    "Codon_position": 2
                }
            ]
        )
        # Strand -, exon, after breakend
        record = VCFRecord(
            "1", 183, "id_21", "A", ["]2:100]A"],
            info={"ANNOT_POS": 183, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 64
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 14
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 20,
                    "Protein_position": 7,
                    "Codon_position": 2
                }
            ]
        )
        # Strand +, exon, after breakend
        record = VCFRecord(
            "1", 183, "id_22", "A", ["]2:100]A"],
            info={"ANNOT_POS": 183, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 64
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 14
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": 7,
                    "Protein_position": 3,
                    "Codon_position": 1
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 20,
                    "Protein_position": 7,
                    "Codon_position": 2
                }
            ]
        )
        # Strand -, exon, before breakend
        record = VCFRecord(
            "1", 183, "id_23", "A", ["A[2:100["],
            info={"ANNOT_POS": 183, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                 {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_01",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_01",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/3&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 64
                },
                {
                    "SYMBOL": "GENE_N01", "Gene": "GENE_I01", "Feature": "TR_02",
                    "Feature_type": "Transcript", "STRAND": "-", "Protein": "PROT_02",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "1/2&5prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None,
                    "CDS_DIST": 14
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "intron",
                    "RNA_ELT_POS": "1/1",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon",
                    "RNA_ELT_POS": "2/3",
                    "CDS_position": 20,
                    "Protein_position": 7,
                    "Codon_position": 2
                }
            ]
        )
        # Strand +, exon, before breakend
        record = VCFRecord(
            "1", 241, "id_24", "A", ["A[2:100["],
            info={"ANNOT_POS": 241, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/2&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "3/3&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )
        # Strand -, exon, after breakend
        record = VCFRecord(
            "1", 241, "id_25", "A", ["]2:100]A"],
            info={"ANNOT_POS": 241, "RNA_FIRST": True, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/2&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "3/3&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )
        # Strand +, exon, after breakend
        record = VCFRecord(
            "1", 241, "id_26", "A", ["]2:100]A"],
            info={"ANNOT_POS": 241, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/2&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "3/3&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )
        # Strand -, exon, before breakend
        record = VCFRecord(
            "1", 241, "id_27", "A", ["A[2:100["],
            info={"ANNOT_POS": 241, "MATEID": "id_02"}
        )
        self.assertEqual(
            sorted(getGeneAnnot(record, genes_by_chr), key=lambda elt: elt["Feature"]),
            [
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_03",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_03",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "2/2&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                },
                {
                    "SYMBOL": "GENE_N02", "Gene": "GENE_I02", "Feature": "TR_04",
                    "Feature_type": "Transcript", "STRAND": "+", "Protein": "PROT_04",
                    "RNA_ELT_TYPE": "exon&utr",
                    "RNA_ELT_POS": "3/3&3prim",
                    "CDS_position": None,
                    "Protein_position": None,
                    "Codon_position": None
                }
            ]
        )

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
