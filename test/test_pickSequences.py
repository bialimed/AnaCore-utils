#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
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
class PickSequences(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_ids = os.path.join(tmp_folder, unique_id + "_in_ids.txt")
        self.tmp_in_seq_1 = os.path.join(tmp_folder, unique_id + "_in.fasta")
        self.tmp_in_seq_2 = os.path.join(tmp_folder, unique_id + "_in.fastq")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out")
        self.tmp_list = [self.tmp_in_ids, self.tmp_in_seq_1, self.tmp_in_seq_2, self.tmp_output]

        # Create input files
        with open(self.tmp_in_ids, "w") as FH_out:
            FH_out.write("""seq1
seq4
seq5
seqMissing""")
        with open(self.tmp_in_seq_1, "w") as FH_out:
            FH_out.write(""">seq1
ATGCATAGTAGTACC
ATAGAGATAGTTTTC
>seq2 desc2
ATGCATAGTAGTACC
ATA
>seq3
ATGCATAGTAGTACC
>seq4

>seq5 desc5
ATGCATAGTAGTACC
ATAGAGATAGTTTTC
GGCAT""")
        with open(self.tmp_in_seq_2, "w") as FH_out:
            FH_out.write("""@seq1
ATGCATAGTAGTACCATAGAGATAGTTTTC
+
##!azaee&&aedd#############aze
@seq2 desc2
ATGCATAGTAGTACCATA
+
##!azaee&&aedd####
@seq3
ATGCATAGTAGTACC
+
##!azaee&&aedd#
@seq4

+

@seq5 desc5
ATGCATAGTAGTACCATAGAGATAGTTTTCGGCAT
+
##!azaee&&aedd#############aze!!:sa""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in self.tmp_list:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testSelectFasta(self):
        # Stringency strict
        cmd = [
            "pickSequences.py",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_1,
            "--output-sequences", self.tmp_output
        ]
        with self.assertRaises(Exception) as context:
            subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [
            "pickSequences.py",
            "--stringency", "lenient",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_1,
            "--output-sequences", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            ">seq1",
            "ATGCATAGTAGTACCATAGAGATAGTTTTC",
            ">seq4",
            "",
            ">seq5 desc5",
            "ATGCATAGTAGTACCATAGAGATAGTTTTCGGCAT"
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )
        # Stringency lenient reprocess to check file handle problems
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )

    def testRemoveFasta(self):
        # Stringency strict
        cmd = [
            "pickSequences.py",
            "--mode", "remove",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_1,
            "--output-sequences", self.tmp_output
        ]
        with self.assertRaises(Exception) as context:
            subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [
            "pickSequences.py",
            "--stringency", "lenient",
            "--mode", "remove",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_1,
            "--output-sequences", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            ">seq2 desc2",
            "ATGCATAGTAGTACCATA",
            ">seq3",
            "ATGCATAGTAGTACC"
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )
        # Stringency lenient reprocess to check file handle problems
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )

    def testSelectFastq(self):
        # Stringency strict
        cmd = [
            "pickSequences.py",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_2,
            "--output-sequences", self.tmp_output
        ]
        with self.assertRaises(Exception) as context:
            subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [
            "pickSequences.py",
            "--stringency", "lenient",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_2,
            "--output-sequences", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            "@seq1",
            "ATGCATAGTAGTACCATAGAGATAGTTTTC",
            "+",
            "##!azaee&&aedd#############aze",
            "@seq4",
            "",
            "+",
            "",
            "@seq5 desc5",
            "ATGCATAGTAGTACCATAGAGATAGTTTTCGGCAT",
            "+",
            "##!azaee&&aedd#############aze!!:sa"
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )
        # Stringency lenient reprocess to check file handle problems
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )

    def testRemoveFastq(self):
        # Stringency strict
        cmd = [
            "pickSequences.py",
            "--mode", "remove",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_2,
            "--output-sequences", self.tmp_output
        ]
        with self.assertRaises(Exception) as context:
            subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Stringency lenient
        cmd = [
            "pickSequences.py",
            "--stringency", "lenient",
            "--mode", "remove",
            "--input-ids", self.tmp_in_ids,
            "--input-sequences", self.tmp_in_seq_2,
            "--output-sequences", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expected = [
            "@seq2 desc2",
            "ATGCATAGTAGTACCATA",
            "+",
            "##!azaee&&aedd####",
            "@seq3",
            "ATGCATAGTAGTACC",
            "+",
            "##!azaee&&aedd#"
        ]
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )
        # Stringency lenient reprocess to check file handle problems
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        observed = None
        with open(self.tmp_output) as FH_obs:
            observed = [elt.strip() for elt in FH_obs.readlines()]
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
