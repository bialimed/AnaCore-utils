#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import pysam
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


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


def bamToSam(in_bam, out_sam):
    with pysam.AlignmentFile(in_bam, "rb") as reader:
        with pysam.AlignmentFile(out_sam, "w", header=reader.header) as writer:
            for record in reader:
                writer.write(record)


class TestSplitBAMByRG(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_sam = os.path.join(tmp_folder, unique_id + "_in.sam")
        self.tmp_in_bam = os.path.join(tmp_folder, unique_id + "_in.bam")
        self.tmp_in_gp = os.path.join(tmp_folder, unique_id + "_in.tsv")
        self.tmp_out_pattern = os.path.join(tmp_folder, unique_id + "_out_{GP}.bam")

        # Create input files
        with open(self.tmp_in_gp, "w") as FH_out:
            FH_out.write("""MGMT	gp1
MLH1	gp2
TEST	gp3""")
        with open(self.tmp_in_sam, "w") as FH_out:
            FH_out.write("""@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:3	LN:198295559
@RG	ID:1	LB:MGMT
@RG	ID:2	LB:MLH1
@RG	ID:3	LB:TEST
M70265:329:000000000-D5GLP:1:1101:15819:1970_Support_ratio:99/100;R1_start:0;R2_start:101	0	3	36993275	42	100M	*	0	0	AGAGTGGATAGTGATTTTTAATGTGTAAGTGTATATTTTTTTAGGTAGTGGGTAGTAGTTGTTTTAGGGAGGGATGAAGAGATTTAGTAATTTATAGAGT	GFFEBGFFFFFFGGGGGGGFFGHHHHHHHHHHHHBGGHHHGGGAFBFAGEEBEGBGDBGGHHHFHGFEFFFFFDHEFFDBABFGHFBGDFGHHHHFF>GF	NM:i:28	MD:Z:4C3C2C3C1C3C1C1C3C1C4C0C2C7C3C5C0C1C2C9C7C0C0C2C2C0C0C1C5	XM:Z:....z...x..z...h.h...z.z.h...z.h....hh..h.......z...x.....xz.h..x.........z.......hhx..h..hhh.x.....	XR:Z:CT	XG:Z:CT	RG:Z:2
M70265:329:000000000-D5GLP:1:1101:17038:2070_Support_ratio:100/100;R1_start:0;R2_start:101	0	3	36993275	42	100M	*	0	0	AGAGTGGATAGTGATTTTTAATGTGTAAGTGTATATTTTTTTAGGTAGTGGGTAGTAGTTGTTTTAGGGAGGGATGAAGAGATTTAGTAATTTATAGAGT	HHHGHHHHHHHHHHGHHHHHHHHGHGHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHGHHFHGHHHGGHHHHHHHHGHHHH	NM:i:28	MD:Z:4C3C2C3C1C3C1C1C3C1C4C0C2C7C3C5C0C1C2C9C7C0C0C2C2C0C0C1C5	XM:Z:....z...x..z...h.h...z.z.h...z.h....hh..h.......z...x.....xz.h..x.........z.......hhx..h..hhh.x.....	XR:Z:CT	XG:Z:CT	RG:Z:2
M70265:329:000000000-D5GLP:1:1101:5862:11858_Support_ratio:100/101;R1_start:0;R2_start:100	0	10	129467210	42	101M	*	0	0	TTGGATATGTTGGGATAGTTTGTGTTTTTAGAATGTTTTGTGTTTTGATGTTTGTAGGTTTTTGTGGTGTGTATTGTTTGTGATTTGGTGAGTGTTTGGGT	F4HHGHHGGHHHHHHGHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHGHHHHHHGHGGHHHHHHHGGG	NM:i:34	MD:Z:0C0C7C5C2C0C0C1C1C0C0C0C5C1C4C2C0C0C2C1C0C0C1C4C0C1C1C4C1C1C0C5C2C11C5	XM:Z:xz.......x.....x..hxz.z.hhhh.....z.h....z..hxz..z.hxz.x....hh.z.z....z.h.xz.....z..h...........x.....	XR:Z:CT	XG:Z:CT	RG:Z:1
M70265:329:000000000-D5GLP:1:1102:19428:6706_Support_ratio:98/101;R1_start:0;R2_start:100	0	10	129467210	42	101M	*	0	0	TTGGATATGTTGGGATAGTTTGTGTTTTTAGAATGTTTTGTGTTTTGATGTTTGTAGGTTTTTGTGGTGTGCATCGTTTGCGATTTGGTGAGTGTTTGGGT	HHHHHHHGGHHHHHHHHHHGHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHG2HHHGHG4G5HG5HHHGG5GGGHGHGHGHHHHHHHHGGG	NM:i:31	MD:Z:0C0C7C5C2C0C0C1C1C0C0C0C5C1C4C2C0C0C2C1C0C0C1C4C0C1C1C4C3C9C11C5	XM:Z:xz.......x.....x..hxz.z.hhhh.....z.h....z..hxz..z.hxz.x....hh.z.z....z.H.xZ.....Z..h...........x.....	XR:Z:CT	XG:Z:CT	RG:Z:1""")

        self.data = [
            {
                "file": self.tmp_out_pattern.replace("{GP}", "gp1"),
                "expected": """@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:3	LN:198295559
@RG	ID:1	LB:MGMT
M70265:329:000000000-D5GLP:1:1101:5862:11858_Support_ratio:100/101;R1_start:0;R2_start:100	0	10	129467210	42	101M	*	0	0	TTGGATATGTTGGGATAGTTTGTGTTTTTAGAATGTTTTGTGTTTTGATGTTTGTAGGTTTTTGTGGTGTGTATTGTTTGTGATTTGGTGAGTGTTTGGGT	F4HHGHHGGHHHHHHGHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHGHHHHHHGHGGHHHHHHHGGG	NM:i:34	MD:Z:0C0C7C5C2C0C0C1C1C0C0C0C5C1C4C2C0C0C2C1C0C0C1C4C0C1C1C4C1C1C0C5C2C11C5	XM:Z:xz.......x.....x..hxz.z.hhhh.....z.h....z..hxz..z.hxz.x....hh.z.z....z.h.xz.....z..h...........x.....	XR:Z:CT	XG:Z:CT	RG:Z:1
M70265:329:000000000-D5GLP:1:1102:19428:6706_Support_ratio:98/101;R1_start:0;R2_start:100	0	10	129467210	42	101M	*	0	0	TTGGATATGTTGGGATAGTTTGTGTTTTTAGAATGTTTTGTGTTTTGATGTTTGTAGGTTTTTGTGGTGTGCATCGTTTGCGATTTGGTGAGTGTTTGGGT	HHHHHHHGGHHHHHHHHHHGHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHG2HHHGHG4G5HG5HHHGG5GGGHGHGHGHHHHHHHHGGG	NM:i:31	MD:Z:0C0C7C5C2C0C0C1C1C0C0C0C5C1C4C2C0C0C2C1C0C0C1C4C0C1C1C4C3C9C11C5	XM:Z:xz.......x.....x..hxz.z.hhhh.....z.h....z..hxz..z.hxz.x....hh.z.z....z.H.xZ.....Z..h...........x.....	XR:Z:CT	XG:Z:CT	RG:Z:1"""
            },
            {
                "file": self.tmp_out_pattern.replace("{GP}", "gp2"),
                "expected": """@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:3	LN:198295559
@RG	ID:2	LB:MLH1
M70265:329:000000000-D5GLP:1:1101:15819:1970_Support_ratio:99/100;R1_start:0;R2_start:101	0	3	36993275	42	100M	*	0	0	AGAGTGGATAGTGATTTTTAATGTGTAAGTGTATATTTTTTTAGGTAGTGGGTAGTAGTTGTTTTAGGGAGGGATGAAGAGATTTAGTAATTTATAGAGT	GFFEBGFFFFFFGGGGGGGFFGHHHHHHHHHHHHBGGHHHGGGAFBFAGEEBEGBGDBGGHHHFHGFEFFFFFDHEFFDBABFGHFBGDFGHHHHFF>GF	NM:i:28	MD:Z:4C3C2C3C1C3C1C1C3C1C4C0C2C7C3C5C0C1C2C9C7C0C0C2C2C0C0C1C5	XM:Z:....z...x..z...h.h...z.z.h...z.h....hh..h.......z...x.....xz.h..x.........z.......hhx..h..hhh.x.....	XR:Z:CT	XG:Z:CT	RG:Z:2
M70265:329:000000000-D5GLP:1:1101:17038:2070_Support_ratio:100/100;R1_start:0;R2_start:101	0	3	36993275	42	100M	*	0	0	AGAGTGGATAGTGATTTTTAATGTGTAAGTGTATATTTTTTTAGGTAGTGGGTAGTAGTTGTTTTAGGGAGGGATGAAGAGATTTAGTAATTTATAGAGT	HHHGHHHHHHHHHHGHHHHHHHHGHGHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHGHHFHGHHHGGHHHHHHHHGHHHH	NM:i:28	MD:Z:4C3C2C3C1C3C1C1C3C1C4C0C2C7C3C5C0C1C2C9C7C0C0C2C2C0C0C1C5	XM:Z:....z...x..z...h.h...z.z.h...z.h....hh..h.......z...x.....xz.h..x.........z.......hhx..h..hhh.x.....	XR:Z:CT	XG:Z:CT	RG:Z:2"""
            },
            {
                "file": self.tmp_out_pattern.replace("{GP}", "gp3"),
                "expected": """@HD	VN:1.0	SO:coordinate
@SQ	SN:1	LN:248956422
@SQ	SN:10	LN:133797422
@SQ	SN:3	LN:198295559
@RG	ID:3	LB:TEST"""
            }
        ]

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_sam, self.tmp_in_gp]:
            if os.path.exists(curr_file):
                os.remove(curr_file)
        for curr_case in self.data:
            if os.path.exists(curr_case["file"]):
                os.remove(curr_case["file"])
            tmp_sam = curr_case["file"].replace(".bam", ".sam")
            if os.path.exists(tmp_sam):
                os.remove(tmp_sam)

    def testSplit(self):
        # Exec
        samToBam(self.tmp_in_sam, self.tmp_in_bam)
        cmd = [
            "splitBAMByRG.py",
            "--input-aln", self.tmp_in_bam,
            "--input-design", self.tmp_in_gp,
            "--output-pattern", self.tmp_out_pattern
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Eval
        for curr_case in self.data:
            tmp_sam = curr_case["file"].replace(".bam", ".sam")
            bamToSam(curr_case["file"], tmp_sam)
            with open(tmp_sam) as reader:
                observed = "".join(reader.readlines())
            self.assertEqual(
                curr_case["expected"].strip(),
                observed.strip()
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
