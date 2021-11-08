#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
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

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
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


class TestSetUMITagFromID(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Data
        self.umi_by_read = {
            "M70265:329:000000000-D5GLP:1:1101:15819:1970": "ATGC+ATTA",
            "M70265:329:000000000-D5GLP:1:1101:17038:2070": "AACC+CTCA",
            "M70265:329:000000000-D5GLP:1:1101:5862:11858": "GGCT+CTAT",
            "M70265:329:000000000-D5GLP:1:1102:19428:6706": "TCGA+TCAC"
        }

        # Temporary files
        self.tmp_in_sam = os.path.join(tmp_folder, unique_id + "_in.sam")
        self.tmp_in_bam = os.path.join(tmp_folder, unique_id + "_in.bam")
        self.tmp_out_bam = os.path.join(tmp_folder, unique_id + "_out.bam")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_bam, self.tmp_in_sam, self.tmp_out_bam]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        for separator in [":", "_"]:
            for umi_tag in ["RX", "MI"]:
                # Prepare SAM
                with open(self.tmp_in_sam, "w") as writer:
                    writer.write("@HD	VN:1.0	SO:coordinate\n")
                    writer.write("@SQ	SN:1	LN:248956422\n")
                    record_idx = 1
                    for read_id, umi_seq in self.umi_by_read.items():
                        writer.write(
                            "{}{}{}	0	1	{}	42	100M	*	0	0	AGAGTGGATAGTGATTTTTAATGTGTAAGTGTATATTTTTTTAGGTAGTGGGTAGTAGTTGTTTTAGGGAGGGATGAAGAGATTTAGTAATTTATAGAGT	GFFEBGFFFFFFGGGGGGGFFGHHHHHHHHHHHHBGGHHHGGGAFBFAGEEBEGBGDBGGHHHFHGFEFFFFFDHEFFDBABFGHFBGDFGHHHHFF>GF	NM:i:28	MD:Z:4C3C2C3C1C3C1C1C3C1C4C0C2C7C3C5C0C1C2C9C7C0C0C2C2C0C0C1C5	XM:Z:....z...x..z...h.h...z.z.h...z.h....hh..h.......z...x.....xz.h..x.........z.......hhx..h..hhh.x.....	XR:Z:CT	XG:Z:CT	RG:Z:2\n".format(
                                read_id, separator, umi_seq, 1000 * record_idx
                            )
                        )
                        record_idx += 1
                samToBam(self.tmp_in_sam, self.tmp_in_bam)
                # Exec script
                cmd = [
                    "setUMITagFromID.py",
                    "--umi-separator", separator,
                    "--umi-tag", umi_tag,
                    "--input-aln", self.tmp_in_bam,
                    "--output-aln", self.tmp_out_bam
                ]
                subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
                # Eval
                observed = {}
                with pysam.AlignmentFile(self.tmp_out_bam, "rb") as reader:
                    for curr_read in reader.fetch(until_eof=True):
                        umi_seq = curr_read.get_tag(umi_tag)
                        observed[curr_read.query_name[:-(len(umi_seq) + 1)]] = curr_read.get_tag(umi_tag)
                self.assertEqual(
                    self.umi_by_read,
                    observed
                )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
