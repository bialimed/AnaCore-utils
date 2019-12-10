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
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class filterVCFBySOR(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "filterVCFBySOR.py",
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {
                "expected": HeaderInfoAttr("expected", "Expected filter tag.", type="String", number="1"),
                "SAR": HeaderInfoAttr("SAR", "Number of reads supporting the alternative allele in reverse strand.", type="Integer", number="1"),
                "SAF": HeaderInfoAttr("SAF", "Number of reads supporting the alternative allele in forward strand.", type="Integer", number="1"),
                "SRR": HeaderInfoAttr("SRR", "Number of reads supporting the reference allele in reverse strand.", type="Integer", number="1"),
                "SRF": HeaderInfoAttr("SRF", "Number of reads supporting the reference allele in forward strand.", type="Integer", number="1"),
            }
            FH_var.writeHeader()
            self.variants = [
                # 0.5 alt, 0.5 ref, low DP, alt no bias, ref no bias
                VCFRecord("artificial_chr1", 10, "sub_01", "G", ["T"], None, None, {
                    "SAR": 5,
                    "SAF": 5,
                    "SRR": 5,
                    "SRF": 5,
                    "expected": "PASS"
                }),
                # 0.05 alt, 0.95 ref, good DP, alt no bias, ref no bias
                VCFRecord("artificial_chr1", 20, "sub_02", "G", ["T"], None, None, {
                    "SAR": 5,
                    "SAF": 5,
                    "SRR": 95,
                    "SRF": 95,
                    "expected": "PASS"
                }),
                # 0.05 alt, 0.95 ref, good DP, alt no bias, ref strand bias
                VCFRecord("artificial_chr1", 30, "sub_03", "G", ["T"], None, None, {
                    "SAR": 5,
                    "SAF": 5,
                    "SRR": 150,
                    "SRF": 30,
                    "expected": "PASS"
                }),
                # 0.05 alt, 0.95 ref, good DP, alt strand bias, ref no bias
                VCFRecord("artificial_chr1", 40, "sub_04", "G", ["T"], None, None, {
                    "SAR": 9,
                    "SAF": 1,
                    "SRR": 95,
                    "SRF": 95,
                    "expected": "strandRatioBias"
                }),
                # 0.05 alt, 0.95 ref, good DP, alt strand bias, ref strand bias => no bias
                VCFRecord("artificial_chr1", 50, "sub_05", "G", ["T"], None, None, {
                    "SAR": 9,
                    "SAF": 1,
                    "SRR": 150,
                    "SRF": 30,
                    "expected": "PASS"
                }),
                # 0.5 alt, 0.5 ref, low DP, alt strand bias, ref no bias
                VCFRecord("artificial_chr1", 60, "sub_06", "G", ["T"], None, None, {
                    "SAR": 9,
                    "SAF": 1,
                    "SRR": 5,
                    "SRF": 5,
                    "expected": "strandRatioBias"
                }),
                # 0.29 alt, 0.71 ref, good DP, alt no bias, ref no bias
                VCFRecord("artificial_chr1", 70, "sub_07", "G", ["T"], None, None, {
                    "SAR": 400,
                    "SAF": 600,
                    "SRR": 1400,
                    "SRF": 1000,
                    "expected": "PASS"
                }),
                # 0.71 alt, 0.29 ref, good DP, alt no bias, ref no bias
                VCFRecord("artificial_chr1", 80, "sub_08", "G", ["T"], None, None, {
                    "SAR": 1400,
                    "SAF": 1000,
                    "SRR": 400,
                    "SRF": 600,
                    "expected": "PASS"
                }),
                # 1.0 alt, 0.0 ref, good DP, alt no bias, ref 0 DP
                VCFRecord("artificial_chr1", 90, "sub_09", "G", ["T"], None, None, {
                    "SAR": 1400,
                    "SAF": 1000,
                    "SRR": 0,
                    "SRF": 0,
                    "expected": "PASS"
                }),
                # 1.0 alt, 0.0 ref, good DP, alt no bias, ref 2 DP
                VCFRecord("artificial_chr1", 100, "sub_10", "G", ["T"], None, None, {
                    "SAR": 1400,
                    "SAF": 1000,
                    "SRR": 0,
                    "SRF": 2,
                    "expected": "PASS"
                }),
                # 1.0 alt, 0.0 ref, limit DP, alt no bias, ref 0 DP
                VCFRecord("artificial_chr1", 110, "sub_11", "G", ["T"], None, None, {
                    "SAR": 90,
                    "SAF": 30,
                    "SRR": 0,
                    "SRF": 0,
                    "expected": "PASS"
                }),
                # 1.0 alt, 0.0 ref, limit DP, alt no bias, ref 2 DP
                VCFRecord("artificial_chr1", 120, "sub_12", "G", ["T"], None, None, {
                    "SAR": 90,
                    "SAF": 30,
                    "SRR": 0,
                    "SRF": 2,
                    "expected": "PASS"
                }),
                # 1.0 alt, 0.0 ref, limit DP, alt strand bias, ref 0 DP
                VCFRecord("artificial_chr1", 130, "sub_13", "G", ["T"], None, None, {
                    "SAR": 90,
                    "SAF": 10,
                    "SRR": 0,
                    "SRF": 0,
                    "expected": "strandRatioBias"
                }),
                # 1.0 alt, 0.0 ref, limit DP, alt strand bias, ref 2 DP
                VCFRecord("artificial_chr1", 140, "sub_14", "G", ["T"], None, None, {
                    "SAR": 90,
                    "SAF": 10,
                    "SRR": 0,
                    "SRF": 2,
                    "expected": "strandRatioBias"
                }),
                # 1.0 alt, 0.0 ref, limit DP, alt strand bias, ref 1 DP
                VCFRecord("artificial_chr1", 150, "sub_15", "G", ["T"], None, None, {
                    "SAR": 90,
                    "SAF": 10,
                    "SRR": 1,
                    "SRF": 0,
                    "expected": "PASS"  # It can be discuss: 2.89
                }),
                # 0.04 alt, 0.96 ref, good DP, alt strand bias, ref no bias
                VCFRecord("artificial_chr1", 160, "sub_16", "G", ["T"], None, None, {
                    "SAR": 15,
                    "SAF": 2,
                    "SRR": 200,
                    "SRF": 200,
                    "expected": "strandRatioBias"
                }),
                # 0.04 alt, 0.96 ref, good DP, alt strand bias, ref no bias
                VCFRecord("artificial_chr1", 170, "sub_17", "G", ["T"], None, None, {
                    "SAR": 13,  # 12 => PASS
                    "SAF": 2,
                    "SRR": 200,
                    "SRF": 200,
                    "expected": "strandRatioBias"
                }),
                # 0.04 alt, 0.96 ref, good DP, alt strand bias, ref strand bias => no bias
                VCFRecord("artificial_chr1", 180, "sub_18", "G", ["T"], None, None, {
                    "SAR": 13,
                    "SAF": 2,
                    "SRR": 350,
                    "SRF": 50,
                    "expected": "PASS"
                }),
                # 0.04 alt, 0.96 ref, good DP, alt strand bias, ref strand bias rev => bias
                VCFRecord("artificial_chr1", 190, "sub_19", "G", ["T"], None, None, {
                    "SAR": 13,
                    "SAF": 2,
                    "SRR": 50,
                    "SRF": 350,
                    "expected": "strandRatioBias"
                }),
                # 0.5 alt, 0.5 ref, low DP, alt strand bias, ref no bias
                VCFRecord("artificial_chr1", 200, "sub_20", "G", ["T"], None, None, {
                    "SAR": 14,
                    "SAF": 2,
                    "SRR": 8,
                    "SRF": 8,
                    "expected": "strandRatioBias"
                }),
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testTag(self):
        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = []
        for record in self.variants:
            for alt in record.alt:
                expected.append(record.id + ":" + record.info["expected"])
        observed = []
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id + ":" + record.filter[0])
        self.assertEqual(
            expected,
            observed
        )

    def testRemove(self):
        # Execute command
        subprocess.check_call(self.cmd + ["--mode", "remove"], stderr=subprocess.DEVNULL)

        # Validate results
        expected = []
        for record in self.variants:
            for alt in record.alt:
                if record.info["expected"] == "PASS":
                    expected.append(record.id)
        observed = []
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
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
