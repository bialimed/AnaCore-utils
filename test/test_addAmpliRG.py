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
import pysam
import tempfile
import unittest
from anacore.bed import getAreas, getAreasByChr, BEDRecord

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)

from addAmpliRG import endOffsetIsUsable, getEndOffset, getOffsetPenalty, getSourceRegion, getStartOffset, selectBestSource


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestAnnotBND(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_aln = os.path.join(tmp_folder, unique_id + "_aln.sam")
        self.tmp_design = os.path.join(tmp_folder, unique_id + "_design.bed")
        with open(self.tmp_design, "w") as writer:
            writer.write("""# Overlap
chr1	49	150	ampl_1	.	.	59	140
chr1	69	170	ampl_2	.	.	79	160
# Inclusion
chr2	95	154	ampl_3	.	.	105	144
chr2	98	149	ampl_4	.	.	108	139
# Near full overlap
chr3	96	155	ampl_5	.	.	106	145
chr3	99	157	ampl_6	.	.	109	147""")
        # ref KRAS_CDS_trimmed: ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAGTACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAAT
        with open(self.tmp_aln, "w") as writer:
            writer.write("""@SQ	SN:chr1	LN:540
@SQ	SN:chr2	LN:540
@SQ	SN:chr3	LN:540
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fasta read_R1.fastq
1_fwd_out_before_offset_45-94	0	chr1	45	60	50M	*	0	0	CAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAAT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_fwd_on_before_offset_46-95	0	chr1	46	60	50M	*	0	0	AAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_fwd_clip_in_before_offset_49-88	0	chr1	49	60	10S40M	*	0	0	AACGTTACTAAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:40	AS:i:40	XS:i:0
1_fwd_on_start_50-99	0	chr1	50	60	50M	*	0	0	GTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGAT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_fwd_on_after_offset_54-103	0	chr1	54	60	50M	*	0	0	CTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_fwd_out_after_offset_55-104	0	chr1	55	60	50M	*	0	0	TTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAAC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_rvs_out_before_offset_96-145	16	chr1	96	60	50M	*	0	0	TGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_rvs_on_start_offset_97-146	16	chr1	97	60	50M	*	0	0	GATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_rvs_on_start_101-150	16	chr1	101	60	50M	*	0	0	CAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_rvs_clip_on_start_offset_110-149	16	chr1	110	60	40M10S	*	0	0	AGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACGCATTACGGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:40	AS:i:40	XS:i:0
1_rvs_on_end_offset_105-154	16	chr1	105	60	50M	*	0	0	AATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
1_rvs_out_after_offset_106-155	16	chr1	106	60	50M	*	0	0	ATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_fwd_out_before_offset_91-140	0	chr2	91	60	50M	*	0	0	GAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_fwd_on_before_offset_92-141	0	chr2	92	60	50M	*	0	0	AATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGAT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_fwd_on_start_96-145	0	chr2	96	60	50M	*	0	0	TGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_fwd_clip_in_after_offset_97-136	0	chr2	97	60	10S40M	*	0	0	TCTGACCGTCGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:40	AS:i:40	XS:i:0
3_fwd_on_after_offset_97-146	0	chr2	97	60	50M	*	0	0	GATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_fwd_on_before_offset_98-147	0	chr2	98	60	50M	*	0	0	ATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_fwd_on_after_offset_endPenalty_98-157	0	chr2	98	60	60M	*	0	0	ATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:60	AS:i:60	XS:i:0
4_fwd_on_before_offset_endPenalty_clipAdapter_98-149	0	chr2	98	60	52M8S	*	0	0	ATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACGAAGAAGG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBAAA@@@??	NM:i:0	MD:Z:52	AS:i:52	XS:i:0
4_fwd_on_start_99-148	0	chr2	99	60	50M	*	0	0	TCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_fwd_on_after_offset_100-149	0	chr2	100	60	50M	*	0	0	CCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAAC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_fwd_on_after_offset_101-150	0	chr2	101	60	50M	*	0	0	CAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_fwd_on_after_offset_103-152	0	chr2	103	60	50M	*	0	0	ACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_rvs_on_before_offset_99-148	16	chr2	99	60	50M	*	0	0	TCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_rvs_on_end_100-149	16	chr2	100	60	50M	*	0	0	CCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAAC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_rvs_on_after_offset_101-150	16	chr2	101	60	50M	*	0	0	CAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_rvs_on_before_offset_102-151	16	chr2	102	60	50M	*	0	0	AACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_rvs_on_before_offset_103-152	16	chr2	103	60	50M	*	0	0	ACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_rvs_on_before_offset_104-153	16	chr2	104	60	50M	*	0	0	CAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_rvs_on_end_105-154	16	chr2	105	60	50M	*	0	0	AATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
3_rvs_after_end_155	16	chr2	106	60	50M	*	0	0	ATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
4_rvs_on_after_offset_endPenalty_clipAdapter_98-150	16	chr2	98	60	10S53M	*	0	0	ACTTTCAGGTATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACC	==>>??@@AABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:53	AS:i:53	XS:i:0
3_rvs_on_before_offset_96-151	16	chr2	96	60	56M	*	0	0	TGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:56	AS:i:56	XS:i:0
5_fwd_96	0	chr3	96	60	50M	*	0	0	TGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
5_fwd_97	0	chr3	97	60	50M	*	0	0	GATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
5_fwd_98	0	chr3	98	60	50M	*	0	0	ATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
6_fwd_99	0	chr3	99	60	50M	*	0	0	TCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
6_fwd_100	0	chr3	100	60	50M	*	0	0	CCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAAC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
6_fwd_101	0	chr3	101	60	50M	*	0	0	CAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACC	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:50	AS:i:50	XS:i:0
6_fwd_endPenalty_99-157	0	chr3	99	60	59M	*	0	0	TCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTCT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:59	AS:i:59	XS:i:0
5_fwd_endPenalty_clipAdapter_99-155_10S	0	chr3	99	60	57M10S	*	0	0	TCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGATGGAGAAACCTGTCTGCCGACCAAG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBAA@@????????	NM:i:0	MD:Z:57	AS:i:57	XS:i:0""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_aln, self.tmp_design]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_getStartOffset(self):
        regions_by_id = {elt.name: elt for elt in getAreas(self.tmp_design)}
        expected = {
            "1_fwd_out_before_offset_45-94": -5,
            "1_fwd_on_before_offset_46-95": -4,
            "1_fwd_clip_in_before_offset_49-88": -1,
            "1_fwd_on_start_50-99": 0,
            "1_fwd_on_after_offset_54-103": 4,
            "1_fwd_out_after_offset_55-104": 5,
            "1_rvs_out_before_offset_96-145": -5,
            "1_rvs_on_start_offset_97-146": -4,
            "1_rvs_clip_on_start_offset_110-149": -1,
            "1_rvs_on_start_101-150": 0,
            "1_rvs_on_end_offset_105-154": 4,
            "1_rvs_out_after_offset_106-155": 5
        }
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                if read.query_name in expected:
                    source_name = "ampl_" + read.query_name.split("_")[0]
                    source = regions_by_id[source_name]
                    self.assertEqual(
                        expected[read.query_name],
                        getStartOffset(read, source)
                    )

    def test_getEndOffset(self):
        regions_by_id = {
            "1_fwd": BEDRecord("chr1", 50, 100),
            "1_rvs": BEDRecord("chr1", 100, 150)
        }
        expected = {
            "1_fwd_out_before_offset_45-94": -6,
            "1_fwd_on_before_offset_46-95": -5,
            "1_fwd_clip_in_before_offset_49-88": -12,
            "1_fwd_on_start_50-99": -1,
            "1_fwd_on_after_offset_54-103": 3,
            "1_fwd_out_after_offset_55-104": 4,
            "1_rvs_out_before_offset_96-145": -4,  # 100 - 96
            "1_rvs_on_start_offset_97-146": -3,  # 100 - 97
            "1_rvs_clip_on_start_offset_110-149": 10,
            "1_rvs_on_start_101-150": 1,
            "1_rvs_on_end_offset_105-154": 5,
            "1_rvs_out_after_offset_106-155": 6
        }
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                if read.query_name in expected:
                    source_name = "_".join(read.query_name.split("_")[:2])
                    source = regions_by_id[source_name]
                    self.assertEqual(
                        expected[read.query_name],
                        getEndOffset(read, source)
                    )

    def getOffsetPenalty(self):
        regions_by_chr = getAreasByChr(self.tmp_design)
        expected = {
            "1_fwd_out_before_offset_45-94": 5 * 1 + 0 * 0.98,
            "1_fwd_on_before_offset_46-95": 4 * 1 + 0 * 0.98,
            "1_fwd_clip_in_before_offset_49-88": 1 * 1 + 0 * 0.98,
            "1_fwd_on_start_50-99": 0 * 1 + 0 * 0.98,
            "1_fwd_on_after_offset_54-103": 4 * 0.99 + 0 * 0.98,
            "1_fwd_out_after_offset_55-104": 5 * 0.99 + 0 * 0.98,
            "1_rvs_out_before_offset_96-145": 5 * 1 + 0 * 0.98,
            "1_rvs_on_start_offset_97-146": 4 * 1 + 0 * 0.98,
            "1_rvs_clip_on_start_offset_110-149": 1 * 1 + 0 * 0.98,
            "1_rvs_on_start_101-150": 0 * 1 + 0 * 0.98,
            "1_rvs_on_end_offset_105-154": 4 * 0.99 + 0 * 0.98,
            "1_rvs_out_after_offset_106-155": 5 * 0.99 + 0 * 0.98,
            "4_rvs_on_after_offset_endPenalty_clipAdapter_98-150": 7 * 0.99 + 1 * 0.98,
            "3_rvs_on_before_offset_96-151": 4 * 0.99 + 1 * 0.98,
            "6_fwd_endPenalty_99-157": 0 * 1 + 0 * 0.98,
            "5_fwd_endPenalty_clipAdapter_99-155_10S": 2 * 0.99 + 0 * 0.98
        }
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                if read.query_name in expected:
                    sources = regions_by_chr[read.reference_name]
                    end_is_usable = endOffsetIsUsable(read, sources[0], sources[1])
                    self.assertEqual(
                        expected[read.query_name],
                        (
                            getOffsetPenalty(read, sources[0], end_is_usable),
                            getOffsetPenalty(read, sources[1], end_is_usable)
                        )
                    )

    def test_endOffsetIsUsable(self):
        regions_by_chr = getAreasByChr(self.tmp_design)
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                sources = regions_by_chr[read.reference_name]
                self.assertEqual(
                    "endPenalty" in read.query_name,
                    endOffsetIsUsable(read, sources[0], sources[1])
                )

    def test_selectBestSource(self):
        regions_by_chr = getAreasByChr(self.tmp_design)
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                if read.query_name not in {"3_rvs_on_before_offset_102-151", "3_rvs_on_before_offset_96-151"}:  # Skip non-realistic values for this test (start mismatch without clipping)
                    expected = "ampl_" + read.query_name.split("_")[0]
                    sources = regions_by_chr[read.reference_name]
                    observed = selectBestSource(read, sources)
                    if observed is not None:
                        observed = observed.name
                    self.assertEqual(
                        expected,
                        str(observed)
                    )
                    observed = selectBestSource(read, sources[::-1])
                    if observed is not None:
                        observed = observed.name
                    self.assertEqual(
                        expected,
                        str(observed)
                    )

    def test_getSourceRegion(self):
        regions_by_chr = getAreasByChr(self.tmp_design)
        with pysam.AlignmentFile(self.tmp_aln, "r") as reader:
            for read in reader.fetch():
                if read.query_name not in {"3_rvs_on_before_offset_102-151", "3_rvs_on_before_offset_96-151"}:  # Skip non-realistic values for this test (start mismatch without clipping)
                    expected = None
                    if "_out_" not in read.query_name:
                        expected = "ampl_" + read.query_name.split("_")[0]
                    sources = regions_by_chr[read.reference_name]
                    observed = getSourceRegion(read, sources, 4)
                    if observed is not None:
                        observed = observed.name
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
