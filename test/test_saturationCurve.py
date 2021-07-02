#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestSaturationCurve(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_counts = os.path.join(tmp_folder, unique_id + ".tsv")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.tsv")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_counts, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testNames(self):
        # Input
        with open(self.tmp_counts, "w") as writer:
            data = [
                ["gene_1", 50, 200],
                ["gene_2", 10, 100],
                ["gene_3", 150, 10],
                ["gene_4", 10, 0],
                ["gene_5", 2, 30],
                ["gene_5", 20, 1],
                ["gene_6", 3, 20]
            ]
            for row in data:
                writer.write("{}\t{}\t{}\n".format(*row))
        # Execute command
        cmd = [
            "saturationCurve.py",
            "--number-steps", "5",
            "--samples-names", "splA", "splB",
            "--input-counts", self.tmp_counts,
            "--output-saturation", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Validate results
        expected = [
            ["sample", "nb_sampled", "nb_unique"],
            ['splA', '49', '6'],
            ['splA', '98', '7'],
            ['splA', '147', '7'],
            ['splA', '196', '7'],
            ['splA', '245', '7'],
            ['splB', '72', '5'],
            ['splB', '144', '6'],
            ['splB', '216', '6'],
            ['splB', '288', '6'],
            ['splB', '360', '6']
        ]
        observed = []
        with open(self.tmp_output) as reader:
            for record in reader:
                observed.append(record.strip().split("\t"))
        self.assertEqual(
            expected,
            observed
        )

    def testNamesInFile(self):
        # Input
        with open(self.tmp_counts, "w") as writer:
            data = [
                ["Gene", "splA", "splB"],
                ["gene_1", 50, 200],
                ["gene_2", 10, 100],
                ["gene_3", 150, 10],
                ["gene_4", 10, 0],
                ["gene_5", 2, 30],
                ["gene_5", 20, 1],
                ["gene_6", 3, 20]
            ]
            for row in data:
                writer.write("{}\t{}\t{}\n".format(*row))
        # Execute command
        cmd = [
            "saturationCurve.py",
            "--number-steps", "5",
            "--input-counts", self.tmp_counts,
            "--output-saturation", self.tmp_output
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Validate results
        expected = [
            ["sample", "nb_sampled", "nb_unique"],
            ['splA', '49', '6'],
            ['splA', '98', '7'],
            ['splA', '147', '7'],
            ['splA', '196', '7'],
            ['splA', '245', '7'],
            ['splB', '72', '5'],
            ['splB', '144', '6'],
            ['splB', '216', '6'],
            ['splB', '288', '6'],
            ['splB', '360', '6']
        ]
        observed = []
        with open(self.tmp_output) as reader:
            for record in reader:
                observed.append(record.strip().split("\t"))
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
