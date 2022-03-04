#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.genomicRegion import CDS, Exon, Gene, Protein, Transcript
from anacore.region import Region, RegionList
import os
import pysam
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from shallowsAnalysis import getTranscriptsAnnot, shallowFromAlignment


########################################################################
#
# FUNCTIONS
#
########################################################################
def samToBam(in_sam, out_bam):
    with pysam.AlignmentFile(in_sam, "r") as reader:
        with pysam.AlignmentFile(out_bam, "wb", header=reader.header) as writer:
            for record in reader:
                writer.write(record)


class DepthAnalysis(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_sam = os.path.join(tmp_folder, unique_id + "_aln.sam")
        self.tmp_bam = os.path.join(tmp_folder, unique_id + "_aln.bam")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sam, self.tmp_bam, self.tmp_bam + ".bai"]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testShallowFromAlignment(self):
        """
        art_chr1:
                10        20        30        40        50        60        70        80        90       100       110       120
        123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|12345678
        ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCCTTGACGATACAGCTAAAT
                   *******.************************************************** ******************************************.*********
                TCGTAAACTTCTGGTAGTTGGAGCTGGTGTTTGCCATAAATAATACTAAATCATTTGAAGA
                                                                              ATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCCTTTACGATACAG
                 ------------------------------------------------------------------------------------------------------------------

        art_chr2:
                10        20        30        40        50        60        70        80        90       100       110       120
        123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|12345678
        ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCCTTGACGATACAGCTAAAT               **********************         **********************
               ********************************************        ***************************************************
               AATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCCATAAATAA        CATTTGAAGATATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCC
                 -------------------------------------------------------------------------------------------

        art_chr3:
                10        20        30        40        50        60        70
        123456789|123456789|123456789|123456789|123456789|123456789|123456789|12
        ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAAAT

                 -------------------------------
        """
        with open(self.tmp_sam, "w") as writer:
            writer.write("""@SQ	SN:art_chr1	LN:128
@SQ	SN:art_chr2	LN:128
@SQ	SN:art_chr3	LN:72
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem ref.fasta reads.fasta
read_1	0	art_chr1	12	60	3S58M	*	0	0	TCGTAAACTTCTGGTAGTTGGAGCTGGTGTTTGCCATAAATAATACTAAATCATTTGAAGA	*	NM:i:1	MD:Z:7G50	AS:i:53	XS:i:0
read_2	0	art_chr1	71	60	52M	*	0	0	ATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCCTTTACGATACAG	*	NM:i:1	MD:Z:42G9	AS:i:47	XS:i:0
read_3	0	art_chr2	8	60	44M8D51M	*	0	0	AATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCCATAAATAACATTTGAAGATATTCACCATTATAGAGAACAAATGCGTAGGCAAGAGTGCC	*	NM:i:8	MD:Z:44^TACTAAAT51	AS:i:81	XS:i:0
""")
        samToBam(self.tmp_sam, self.tmp_bam)
        pysam.index(self.tmp_bam)

        class FakeLogger:
            def info(self, msg):
                pass

        selected_regions = RegionList([
            Region(10, 123, None, "art_chr1"),
            Region(10, 100, None, "art_chr2"),
            Region(10, 40, None, "art_chr3"),
        ])

        expected = [
            "art_chr1:10-11",
            "art_chr1:70-70",
            "art_chr1:123-123",
            "art_chr3:10-40",
        ]
        observed = [str(elt) for elt in shallowFromAlignment(self.tmp_bam, selected_regions, "reads", 10, 1, FakeLogger())]
        self.assertEqual(sorted(expected), sorted(observed))

    def testGetTranscriptsAnnot_withUTR_threeExons(self):
        exon_1 = Exon(10, 40, "+", "chr1", "fwd_exon_1")
        exon_2 = Exon(91, 150, "+", "chr1", "fwd_exon_2")
        exon_3 = Exon(201, 361, "+", "chr1", "fwd_exon_3")
        cds_2 = CDS(110, 150, "+", "chr1", "fwd_cds_2")
        cds_3 = CDS(201, 246, "+", "chr1", "fwd_cds_3")
        gene_1 = Gene(10, 350, None, "chr1", "gene_1", {"id": "g_1"})
        transcrit_1 = Transcript(None, None, None, "chr1", "transcrit_1", {"id": "tr_1"}, parent=gene_1, children=[exon_1, exon_2, exon_3])
        protein_1 = Protein(None, None, None, "chr1", "protein_1", children=[cds_2, cds_3], transcript=transcrit_1)
        queries = [
            Region(80, 100, None, "chr1", "query_1", {"desc": "starts before exon_2 ; ends in exon_2 before CDS."}),
            Region(80, 115, None, "chr1", "query_2", {"desc": "starts before exon_2 ; ends in exon_2 in CDS."}),
            Region(94, 110, None, "chr1", "query_3", {"desc": "starts in exon_2 before CDS ; ends in exon_2 at start of CDS."}),
            Region(100, 180, None, "chr1", "query_4", {"desc": "starts in exon_2 before CDS ; ends after exon_2."}),
            Region(91, 150, None, "chr1", "query_5", {"desc": "starts in exon_2 before CDS ; ends at end of exon_2 in CDS."}),
            Region(80, 170, None, "chr1", "query_6", {"desc": "starts before exon_2 ; ends after exon_2."}),
            Region(80, 230, None, "chr1", "query_7", {"desc": "starts before exon_2 ; ends in exon_3 in CDS."}),
            Region(100, 400, None, "chr1", "query_8", {"desc": "starts in exon_2 before CDS ; ends after exon_3."}),
            Region(100, 250, None, "chr1", "query_9", {"desc": "starts in exon_2 before CDS ; ends in exon_3 after CDS."}),
            Region(80, 370, None, "chr1", "query_10", {"desc": "starts before exon_2 ; ends after exon_3."}),
            Region(143, 230, None, "chr1", "query_11", {"desc": "starts in exon_2 in CDS ; ends in exon_3 in CDS."}),
            Region(110, 246, None, "chr1", "query_12", {"desc": "starts in exon_2 at start of CDS ; ends in exon_3 at end of CDS."})
        ]

        # Expected forward 3 exons
        expected = {
            "query_1": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": None,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": None
            },
            "query_2": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 2
            },
            "query_3": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 1
            },
            "query_4": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 14
            },
            "query_5": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 14
            },
            "query_6": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 1,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 14
            },
            "query_7": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 1,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 24
            },
            "query_8": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_9": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_10": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 1,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_11": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 12,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 24
            },
            "query_12": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "1"}.items():
                query_res[key] = val
        # Apply forward strand
        for exon in transcrit_1.children:
            exon.strand = "+"
        for cds in protein_1.children:
            cds.strand = "+"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Assert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

        # Expected reverse 3 exons
        expected = {
            "query_1": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": None,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": None
            },
            "query_2": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 28,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 29
            },
            "query_3": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 29,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_4": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 16,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_5": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 16,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_6": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 16,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 29
            },
            "query_7": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 6,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 29
            },
            "query_8": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_9": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
            "query_10": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 29
            },
            "query_11": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 6,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_12": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 29
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "-1"}.items():
                query_res[key] = val
        # Apply reverse strand
        for exon in transcrit_1.children:
            exon.strand = "-"
        for cds in protein_1.children:
            cds.strand = "-"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Assert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

    def testGetTranscriptsAnnot_withoutUTR_threeExons(self):
        exon_1 = Exon(10, 40, "+", "chr1", "fwd_exon_1")
        exon_2 = Exon(91, 150, "+", "chr1", "fwd_exon_2")
        exon_3 = Exon(201, 361, "+", "chr1", "fwd_exon_3")
        cds_1 = CDS(10, 40, "+", "chr1", "fwd_cds_1")
        cds_2 = CDS(91, 150, "+", "chr1", "fwd_cds_2")
        cds_3 = CDS(201, 361, "+", "chr1", "fwd_cds_3")
        gene_1 = Gene(10, 350, None, "chr1", "gene_1", {"id": "g_1"})
        transcrit_1 = Transcript(None, None, None, "chr1", "transcrit_1", {"id": "tr_1"}, parent=gene_1, children=[exon_1, exon_2, exon_3])
        protein_1 = Protein(None, None, None, "chr1", "protein_1", children=[cds_1, cds_2, cds_3], transcript=transcrit_1)
        queries = [
            Region(80, 100, None, "chr1", "query_1", {"desc": "starts before exon_2 ; ends in exon_2."}),
            Region(100, 180, None, "chr1", "query_2", {"desc": "starts in exon_2 ; ends after exon_2."}),
            Region(91, 150, None, "chr1", "query_3", {"desc": "starts at the start of exon_2 ; ends at the end of exon_2."}),
            Region(80, 170, None, "chr1", "query_4", {"desc": "starts before exon_2 ; ends after exon_2."}),
            Region(80, 230, None, "chr1", "query_5", {"desc": "starts before exon_2 ; ends in exon_3."}),
            Region(100, 400, None, "chr1", "query_6", {"desc": "starts in exon_2 ; ends after exon_3."}),
            Region(100, 250, None, "chr1", "query_7", {"desc": "starts in exon_2 ; ends in exon_3."}),
            Region(80, 370, None, "chr1", "query_8", {"desc": "starts before exon_2 ; ends after exon_3."}),
            Region(90, 151, None, "chr1", "query_9", {"desc": "starts just before exon_2 ; ends just after exon_2."})
        ]

        # Expected forward 3 exons
        expected = {
            "query_1": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 11,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 14
            },
            "query_2": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 14,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 31
            },
            "query_3": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 11,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 31
            },
            "query_4": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 11,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 31
            },
            "query_5": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 11,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 41
            },
            "query_6": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 14,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 84
            },
            "query_7": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 14,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 47
            },
            "query_8": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 11,
                "end_EXON": "3/3",
                "end_INTRON": None,
                "end_Protein_position": 84
            },
            "query_9": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 11,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 31
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "1"}.items():
                query_res[key] = val
        # Apply forward strand
        for exon in transcrit_1.children:
            exon.strand = "+"
        for cds in protein_1.children:
            cds.strand = "+"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Assert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

        # Expected reverse 3 exons
        expected = {
            "query_1": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 71,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 74
            },
            "query_2": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 54,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 71
            },
            "query_3": {
                "start_EXON": "2/3",
                "start_INTRON": None,
                "start_Protein_position": 54,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 74
            },
            "query_4": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 54,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 74
            },
            "query_5": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 44,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 74
            },
            "query_6": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 71
            },
            "query_7": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 38,
                "end_EXON": "2/3",
                "end_INTRON": None,
                "end_Protein_position": 71
            },
            "query_8": {
                "start_EXON": "1/3",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 74
            },
            "query_9": {
                "start_EXON": None,
                "start_INTRON": "1/2",
                "start_Protein_position": 54,
                "end_EXON": None,
                "end_INTRON": "2/2",
                "end_Protein_position": 74
            }
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "-1"}.items():
                query_res[key] = val
        # Apply reverse strand
        for exon in transcrit_1.children:
            exon.strand = "-"
        for cds in protein_1.children:
            cds.strand = "-"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Asert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

    def testGetTranscriptsAnnot_withUTR_oneExon(self):
        exon_1 = Exon(91, 150, "+", "chr1", "exon_2")
        cds_1 = CDS(94, 147, "+", "chr1", "cds_1")
        gene_1 = Gene(10, 350, None, "chr1", "gene_1", {"id": "g_1"})
        transcrit_1 = Transcript(None, None, None, "chr1", "transcrit_1", {"id": "tr_1"}, parent=gene_1, children=[exon_1])
        protein_1 = Protein(None, None, None, "chr1", "protein_2", children=[cds_1], transcript=transcrit_1)
        queries = [
            Region(80, 160, None, "chr1", "query_1", {"desc": "starts before exon_1 ; ends after exon_1."}),
            Region(91, 150, None, "chr1", "query_2", {"desc": "starts at start of exon_1 before CDS ; ends at end of exon_1 after CDS."}),
            Region(94, 147, None, "chr1", "query_3", {"desc": "starts in exon_1 at start of CDS ; ends in exon_1 at end of CDS."}),
            Region(92, 148, None, "chr1", "query_4", {"desc": "starts in exon_1 before CDS ; ends in exon_1 after CDS."}),
            Region(100, 110, None, "chr1", "query_5", {"desc": "starts in exon_1 in CDS ; ends in exon_1 in CDS."}),
            Region(80, 100, None, "chr1", "query_6", {"desc": "starts before exon_1 ; ends in exon_1 in CDS."}),
            Region(92, 110, None, "chr1", "query_7", {"desc": "starts in exon_1 before CDS ; ends in exon_1 in CDS."}),
            Region(110, 200, None, "chr1", "query_8", {"desc": "starts in exon_1 in CDS ; ends after exon_1."}),
            Region(110, 148, None, "chr1", "query_9", {"desc": "starts in exon_1 in CDS ; ends in exon_1 after CDS."}),
        ]

        # Expected forward 1 exon
        expected = {
            "query_1": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_2": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_3": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_4": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_5": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 3,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 6
            },
            "query_6": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 3
            },
            "query_7": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 6
            },
            "query_8": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 6,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_9": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 6,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "1"}.items():
                query_res[key] = val
        # Apply forward strand
        for exon in transcrit_1.children:
            exon.strand = "+"
        for cds in protein_1.children:
            cds.strand = "+"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Asert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

        # Expected reverse 1 exon
        expected = {
            "query_1": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_2": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_3": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_4": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_5": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 13,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 16
            },
            "query_6": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 16,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_7": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 13,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 18
            },
            "query_8": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 13
            },
            "query_9": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 13
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "-1"}.items():
                query_res[key] = val
        # Apply reverse strand
        for exon in transcrit_1.children:
            exon.strand = "-"
        for cds in protein_1.children:
            cds.strand = "-"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Asert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

    def testGetTranscriptsAnnot_withoutUTR_oneExon(self):
        exon_1 = Exon(91, 150, "+", "chr1", "exon_2")
        cds_1 = CDS(91, 150, "+", "chr1", "cds_1")
        gene_1 = Gene(10, 350, None, "chr1", "gene_1", {"id": "g_1"})
        transcrit_1 = Transcript(None, None, None, "chr1", "transcrit_1", {"id": "tr_1"}, parent=gene_1, children=[exon_1])
        protein_1 = Protein(None, None, None, "chr1", "protein_2", children=[cds_1], transcript=transcrit_1)
        queries = [
            Region(80, 160, None, "chr1", "query_1", {"desc": "starts before exon_1 ; ends after exon_1."}),
            Region(91, 150, None, "chr1", "query_2", {"desc": "starts at start of exon_1 ; ends at end of exon_1."}),
            Region(100, 110, None, "chr1", "query_3", {"desc": "starts in exon_1 ; ends in exon_1."}),
            Region(80, 100, None, "chr1", "query_4", {"desc": "starts before exon_1 ; ends in exon_1."}),
            Region(110, 200, None, "chr1", "query_5", {"desc": "starts in exon_1 ; ends after exon_1."}),
        ]

        # Expected forward 1 exon
        expected = {
            "query_1": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
            "query_2": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
            "query_3": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 4,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 7
            },
            "query_4": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 4
            },
            "query_5": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 7,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "1"}.items():
                query_res[key] = val
        # Apply forward strand
        for exon in transcrit_1.children:
            exon.strand = "+"
        for cds in protein_1.children:
            cds.strand = "+"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Asert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)

        # Expected reverse 1 exon
        expected = {
            "query_1": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
            "query_2": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
            "query_3": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 14,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 17
            },
            "query_4": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 17,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 20
            },
            "query_5": {
                "start_EXON": "1/1",
                "start_INTRON": None,
                "start_Protein_position": 1,
                "end_EXON": "1/1",
                "end_INTRON": None,
                "end_Protein_position": 14
            },
        }
        for query_name, query_res in expected.items():
            for key, val in {"SYMBOL": "gene_1", "Gene": "g_1", "Feature": "tr_1", "Feature_type": "Transcript", "STRAND": "-1"}.items():
                query_res[key] = val
        # Apply reverse strand
        for exon in transcrit_1.children:
            exon.strand = "-"
        for cds in protein_1.children:
            cds.strand = "-"
        transcrit_1.sortChildren()
        protein_1.sortChildren()
        # Asert
        for curr_query in queries:
            annotations = getTranscriptsAnnot(curr_query, [transcrit_1])
            self.assertEqual([expected[curr_query.name]], annotations)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
