#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

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

from depthsPanel import depthsTargets, getTargets


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


class DepthsPanel(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_tgt = os.path.join(tmp_folder, unique_id + "_tgt.tsv")
        self.tmp_sam = os.path.join(tmp_folder, unique_id + "_aln.sam")
        self.tmp_bam = os.path.join(tmp_folder, unique_id + "_aln.bam")
        # Write targets
        with open(self.tmp_tgt, "w") as writer:
            writer.write("""name	locations	transcript	description
KRAS	chr1:5-12,chr1:20-36	NM_0505	Exon 1 and 2
TERT	chr2:1-6	NM_77878	Promoter
spl_tracking	chr1:10-10,chr2:8-8 		SNP for identity tracking
missing	chr3:9-13 		""")
        # Write alignments
        with open(self.tmp_sam, "w") as writer:
            writer.write("""@SQ	SN:chr1	LN:60
@SQ	SN:chr2	LN:40
@SQ	SN:chr3	LN:28
@PG	ID:bwa	PN:bwa VN:0.7.10-r789 CL:bwa mem ref.fasta reads.fasta
read1	0	chr1	1	60	34M	*	0	0	ATGACTGAATATAAACTTGTGGTAGTTGGAGCTG	*	NM:i:0	MD:Z:34	AS:i:34	XS:i:0
read3	0	chr1	1	60	59M	*	0	0	ATGACTGAATATAAACTTGTGGTAATTGGAGCTGGTGTTTGCCATAAATAATACTAAAT	*	NM:i:1	MD:Z:24G34	AS:i:54	XS:i:0
read2	0	chr1	9	60	35M	*	0	0	ATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCC	*	NM:i:0	MD:Z:35	AS:i:35	XS:i:0
read4	0	chr2	1	60	31M	*	0	0	ATTTGAAGATATTCACCATTATAGAGAACAA	*	NM:i:0	MD:Z:31	AS:i:31	XS:i:0
""")
        # chr1:
        #          10        20        30        40        50
        # 123456789|123456789|123456789|123456789|123456789|123456789|
        # ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCCATAAATAATACTAAATC
        #     ........       .................
        #
        # Aln:
        # ATGACTGAATATAAACTTGTGGTAGTTGGAGCTG
        #         ATATAAACTTGTGGTAGTTGGAGCTGGTGTTTGCC
        # ATGACTGAATATAAACTTGTGGTAATTGGAGCTGGTGTTTGCCATAAATAATACTAAAT  # mistmacth in second target
        #
        # chr2:
        #          10        20        30
        # 123456789|123456789|123456789|123456789|
        # ATTTGAAGATATTCACCATTATAGAGAACAAATGCGTAGG
        # ......
        #
        # Aln:
        # ATTTGAAGATATTCACCATTATAGAGAACAA
        #
        # chr3:
        #          10        20
        # 123456789|123456789|12345678
        # CAAGAGTGCCTTGACGATACAGCTAAAT

        samToBam(self.tmp_sam, self.tmp_bam)
        pysam.index(self.tmp_bam)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_tgt, self.tmp_sam, self.tmp_bam, self.tmp_bam + ".bai"]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_depthsTargets(self):
        class FakeLogger:
            def info(self, msg):
                pass
        expected = [
            {1: 0, 2: 0, 3: 6},
            {1: 0, 2: 6, 3: 6},
            {1: 0, 2: 1, 3: 1},
            {1: 5, 2: 5, 3: 5}
        ]
        targets = getTargets(self.tmp_tgt)
        depthsTargets(self.tmp_bam, targets, "reads", [1, 2, 3], FakeLogger())
        observed = [elt["under_threshold"] for elt in targets]
        self.assertEqual(expected, observed)

    def test_getTargets(self):
        expected = [
            {
                "name": "KRAS",
                "locations": ["chr1:5-12", "chr1:20-36"],
                "size": 8 + 17,
                "metadata": {"transcript": "NM_0505", "description": "Exon 1 and 2"}
            },
            {
                "name": "TERT",
                "locations": ["chr2:1-6"],
                "size": 6,
                "metadata": {"transcript": "NM_77878", "description": "Promoter"}
            },
            {
                "name": "spl_tracking",
                "locations": ["chr1:10-10", "chr2:8-8"],
                "size": 2,
                "metadata": {"description": "SNP for identity tracking"}
            },
            {
                "name": "missing",
                "locations": ["chr3:9-13"],
                "size": 5,
                "metadata": {}
            }
        ]
        observed = getTargets(self.tmp_tgt)
        for curr in observed:
            curr["locations"] = [str(elt) for elt in curr["locations"]]
        self.assertEqual(expected, observed)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
