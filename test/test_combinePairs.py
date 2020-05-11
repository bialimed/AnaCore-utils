#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
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
class CombinePairs(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_R1 = os.path.join(tmp_folder, unique_id + "_R1.fastq.gz")
        self.tmp_R2 = os.path.join(tmp_folder, unique_id + "_R2.fastq.gz")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_combined.fastq")

        # Exec command
        self.cmd = [
            "combinePairs.py",
            "--input-R1", self.tmp_R1,
            "--input-R2", self.tmp_R2,
            "--output-combined", self.tmp_output
        ]

        # Create input file R1
        a_content = """@seq_1	Empty

+

@seq_2	Too short
CACTAAATAGCAGCATTCA
+
HGFEDCBA@?>=<;:9876
@seq_3	Fragment length > reads length and errors 96:A>C,116:A>G
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCACATGCCAAAGTAGCTGATGTGTTGGACGTTCTAAATGAGGTAGATGAATATTCTG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC#BBBAA@@@@@@@@@=<<<<;<<<<;;;;;;;;;;;;;;;;::::::::::::98
@seq_4	Fragment length < reads length and errors 80:T>C
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGACGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTAGCGAATTTCGACGATCGTTGCATTAACTGTAGTA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCBBBBAA@@@@@@@@@<<<<<<<<<<;;;;;;;;;;;;;;;;::::::::::::98
@seq_5	Too short overlap: on poly-A size 15
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACAAAAAAAAAAAAAAAA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCC
@seq_6	overlap on last 20nt and errors 80:A>T
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACAAAATAAAAAAAAAAATTGT
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCBBBB"""
        with open(self.tmp_R1, "w") as FH_out:
            FH_out.write(a_content)
        # @seq_4	Fragment length < reads length
        # Create input file R2
        b_content = """@seq_1	Empty

+

@seq_2	Too short
TGAATGCTGCTATTTAGTG
+
FEDCBA@?>=<;:987654
@seq_3	Fragment length > reads length and errors 87:A>T,149:A>G,150:C>G
TTTTACATATTAAAGCCTCATGAGGATCACTGGCCAGTAAGTCTATTTTCTCTGAAGAACCAGAATATTCATCTACCTCATTTAGATCGTCCAATACATCAGCTACTTTGGCATTTGATTCAGACTCCCCATCATGTGAGTCATCAGAGG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=CCCCCCCCBBBBAA@@@@@@@@@=<<<<<<<<<;;;;;;;;;;;;;;;;::::::::::::98
@seq_4	Fragment length < reads length and errors 73:T>C,110:T>A
TACATCAGCTACTTTGGCATTTGATTCAGACTCCCCATCATGTGAGTCATCAGAACCTAACAGTTCATCACTCCTGGAAAACCACTCATTAACTTTCTGAATGCTGCTAATTAGTGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCBBBBAA@@@@@@@@9<<<<<<<<<<;;;;;;;;;;;;;;;;::::::::::::98
@seq_5	Too short overlap: on poly-A size 15
CCTCATTTAGAACGTCCAATACATCAGCTACTTTGGCATTTGATTCAGACTCCCCATCATGTGAGTCATCAGAACTTTTTTTTTTTTTTTT
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCC
@seq_6	overlap on last 20nt
CCTCATTTAGAACGTCCAATACATCAGCTACTTTGGCATTTGATTCAGACTCCCCATCATGTGAGTCATCAGAACACAATTTTTTTTTTTTTTTT
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCBB@@@"""
        with open(self.tmp_R2, "w") as FH_out:
            FH_out.write(b_content)


    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_R1, self.tmp_R2, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """@seq_3 Support_ratio:85/90;R1_start:60;R2_start:0
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGATCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCTCATGAGGCTTTAATATGTAAAA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC<BBBAA@@@@@@@@AABBBBCCCCCCCC=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH
@seq_4 Support_ratio:113/116;R1_start:0;R2_start:34
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH
@seq_6 Support_ratio:19/20;R1_start:75;R2_start:0
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACAAAAAAAAAAAAAAAATTGTGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGTTCTAAATGAGG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH"""
        observed = ""
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            sorted(expected.split("\n")),
            sorted(observed.split("\n"))
        )


    def testMaxFragLength(self):
        # Execute command
        custom_cmd = [elt for elt in self.cmd]
        custom_cmd.extend(["--max-frag-length", "175"])
        subprocess.check_call(custom_cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """@seq_4 Support_ratio:113/116;R1_start:0;R2_start:34
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH
@seq_6 Support_ratio:19/20;R1_start:75;R2_start:0
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACAAAAAAAAAAAAAAAATTGTGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGTTCTAAATGAGG
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBCCCCCCCCCCCCCCCCCACCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH"""
        observed = ""
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            sorted(expected.split("\n")),
            sorted(observed.split("\n"))
        )


    def testMinFragLength(self):
        # Execute command
        custom_cmd = [elt for elt in self.cmd]
        custom_cmd.extend(["--min-frag-length", "180"])
        subprocess.check_call(custom_cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """@seq_3 Support_ratio:85/90;R1_start:60;R2_start:0
CACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTGATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATTGGACGATCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCTCATGAGGCTTTAATATGTAAAA
+
HGHHHHGGGGGGGGGGGGGGGGFFFFFFFGFFFEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC<BBBAA@@@@@@@@AABBBBCCCCCCCC=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEFFFGFFFFFFFGGGGGGGGGGGGGGGGHHHHGH"""
        observed = ""
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        self.assertEqual(
            sorted(expected.split("\n")),
            sorted(observed.split("\n"))
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
