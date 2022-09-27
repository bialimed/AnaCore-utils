#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import subprocess
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestStatusToAnnot(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_annot.tsv")
        self.tmp_status = os.path.join(tmp_folder, unique_id + "_status.tsv")
        with open(self.tmp_status, "w") as writer:
            writer.write("""sample	BAT25	BAT26	NR27	status
splA	MSS	MSS	MSS	MSS
splB	MSI	MSI	MSI	MSI
splC	MSI	MSS	MSI	MSI""")
        self.tmp_status_by_locus = os.path.join(tmp_folder, unique_id + "_statusByLocus.tsv")
        with open(self.tmp_status_by_locus, "w") as writer:
            writer.write("""sample	4:54732007-54732108	2:47414383-47414484	11:102322740-102322841	status
splA	MSS	MSS	MSS	MSS
splB	MSI	MSI	MSI	MSI
splC	MSI	MSS	MSI	MSI""")
        self.tmp_targets = os.path.join(tmp_folder, unique_id + "_targets.tsv")
        with open(self.tmp_targets, "w") as writer:
            writer.write("""2	47414383	47414484	BAT26
4	54732007	54732108	BAT25
11	102322740	102322841	NR27""")
        self.expect = """sample	locus_position	method_id	key	value	type
splA	2:47414383-47414484	model	status	MSS	str
splA	4:54732007-54732108	model	status	MSS	str
splA	11:102322740-102322841	model	status	MSS	str
splB	2:47414383-47414484	model	status	MSI	str
splB	4:54732007-54732108	model	status	MSI	str
splB	11:102322740-102322841	model	status	MSI	str
splC	2:47414383-47414484	model	status	MSS	str
splC	4:54732007-54732108	model	status	MSI	str
splC	11:102322740-102322841	model	status	MSI	str""".split("\n")

    def tearDown(self):
        for curr_file in [self.tmp_out, self.tmp_status, self.tmp_status_by_locus, self.tmp_targets]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_by_locus(self):
        cmd = [
            "microsatStatusToAnnot.py",
            "--locus-id",
            "--input-status", self.tmp_status_by_locus,
            "--input-targets", self.tmp_targets,
            "--output-annotations", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            obs = [elt.strip() for elt in reader.readlines()]
        self.assertEqual(sorted(obs), sorted(self.expect))

    def test_by_name(self):
        cmd = [
            "microsatStatusToAnnot.py",
            "--input-status", self.tmp_status,
            "--input-targets", self.tmp_targets,
            "--output-annotations", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            obs = [elt.strip() for elt in reader.readlines()]
        self.assertEqual(sorted(obs), sorted(self.expect))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
