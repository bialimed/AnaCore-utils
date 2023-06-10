#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess
from anacore.vcf import VCFRecord, HeaderInfoAttr
from anacore.annotVcf import AnnotVCFIO

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestFilterAnnotVCF(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_var_filters = os.path.join(tmp_folder, unique_id + "_varFilters.json")
        self.tmp_annot_filters = os.path.join(tmp_folder, unique_id + "_annFilters.json")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Command
        self.cmd = [
            "filterAnnotVCF.py",
            "--input-variants", self.tmp_variants,
            "--output-variants", self.tmp_output
        ]

        # Create filters
        with open(self.tmp_var_filters, "w") as FH_filter:
            FH_filter.write("""{
    "class": "FiltersCombiner",
    "operator": "or",
    "filters": [
        {
            "class": "Filter",
            "getter": "filter",
            "action": "select",
            "aggregator": "ratio:1",
            "operator": "!=",
            "values": "CSQ"
        }, {
            "class": "Filter",
            "getter": "chrom",
            "action": "select",
            "operator": "==",
            "values": "artificial_chr2"
        }
    ]
}""")
        with open(self.tmp_annot_filters, "w") as FH_filter:
            FH_filter.write("""{
    "class": "Filter",
    "getter": "FILTER",
    "action": "select",
    "operator": "==",
    "values": "PASS"
}""")

        # Create VCF
        with AnnotVCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.ANN_titles = ["Allele", "id", "is_filtered", "FILTER"]
            FH_var.info = {
                "ANN": HeaderInfoAttr("ANN", "Consequence annotations from Ensembl VEP. Format: Allele|id|is_filtered|FILTER.", type="String", number="."),
                "is_filtered": HeaderInfoAttr("is_filtered", "The expected result.", type="Integer", number="1")
            }
            FH_var.writeHeader()
            self.variants = [
                VCFRecord(
                    "artificial_chr1", 10, "alt_00", "G", ["T"], None, ["PASS"],
                    {
                        "is_filtered": 0
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_01", "G", ["T"], None, ["CSQ"],
                    {
                        "is_filtered": 1
                    }
                ),
                VCFRecord(
                    "artificial_chr2", 10, "alt_02", "G", ["T"], None, ["CSQ"],
                    {
                        "is_filtered": 0,  # Proctected
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_03", "G", ["T"], None, ["PASS"],
                    {
                        "ANN": [{"Allele": "T", "id": "ann_00", "FILTER": "PASS", "is_filtered": 0}],
                        "is_filtered": 0
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_04", "G", ["T"], None, ["PASS"],
                    {
                        "ANN": [{"Allele": "C", "id": "ann_01", "FILTER": "ANN.COLLOC", "is_filtered": 1}],
                        "is_filtered": 0
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_05", "G", ["T"], None, ["CSQ"],
                    {
                        "ANN": [{"Allele": "C", "id": "ann_02", "FILTER": "ANN.COLLOC", "is_filtered": 1}],
                        "is_filtered": 1
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_06", "G", ["T"], None, ["CSQ"],
                    {
                        "ANN": [{"Allele": "T", "id": "ann_03", "FILTER": "PASS", "is_filtered": 0}],
                        "is_filtered": 1
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_07", "G", ["T"], None, ["PASS"],
                    {
                        "ANN": [
                            {"Allele": "T", "id": "ann_04", "FILTER": "PASS", "is_filtered": 0},
                            {"Allele": "C", "id": "ann_05", "FILTER": "ANN.COLLOC", "is_filtered": 1},
                        ],
                        "is_filtered": 0
                    }
                ),
                VCFRecord(
                    "artificial_chr1", 10, "alt_08", "G", ["T"], None, ["PASS"],
                    {
                        "ANN": [
                            {"Allele": "T", "id": "ann_06", "FILTER": "ANN.popAF", "is_filtered": 1},
                            {"Allele": "C", "id": "ann_07", "FILTER": "ANN.COLLOC&ANN.popAF", "is_filtered": 1},
                        ],
                        "is_filtered": 0
                    }
                ),
                VCFRecord(
                    "artificial_chr2", 10, "alt_09", "G", ["T"], None, ["CSQ"],
                    {
                        "ANN": [
                            {"Allele": "T", "id": "ann_08", "FILTER": "ANN.popAF", "is_filtered": 1},
                            {"Allele": "C", "id": "ann_09", "FILTER": "ANN.COLLOC&ANN.popAF", "is_filtered": 1},
                        ],
                        "is_filtered": 0  # Protected
                    }
                ),
                VCFRecord(
                    "artificial_chr2", 10, "alt_10", "G", ["T"], None, ["CSQ"],
                    {
                        "ANN": [
                            {"Allele": "T", "id": "ann_10", "FILTER": "PASS", "is_filtered": 0},
                            {"Allele": "C", "id": "ann_11", "FILTER": "ANN.COLLOC&ANN.popAF", "is_filtered": 1},
                        ],
                        "is_filtered": 0  # Protected
                    }
                )
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_var_filters, self.tmp_annot_filters, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testMissingFilterException(self):
        with self.assertRaises(subprocess.CalledProcessError) as context:
            subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

    def testResultsFilterVar(self):
        cmd = self.cmd + ["--input-filters-variants", self.tmp_var_filters]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            if record.info["is_filtered"] == 0:
                expected.append(record.id)
                if "ANN" in record.info:
                    for curr_ann in record.info["ANN"]:
                        expected.append(curr_ann["id"])
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
                if "ANN" in record.info:
                    for curr_ann in record.info["ANN"]:
                        observed.append(curr_ann["id"])
        self.assertEqual(observed, expected)

    def testResultsFilterAnn(self):
        cmd = self.cmd + ["--input-filters-annotations", self.tmp_annot_filters]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            expected.append(record.id)
            if "ANN" in record.info:
                for curr_ann in record.info["ANN"]:
                    if curr_ann["is_filtered"] == 0:
                        expected.append(curr_ann["id"])
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
                if "ANN" in record.info:
                    for curr_ann in record.info["ANN"]:
                        observed.append(curr_ann["id"])

        self.assertEqual(observed, expected)

    def testResultsFilterAnnAndVar(self):
        cmd = self.cmd + ["--input-filters-variants", self.tmp_var_filters, "--input-filters-annotations", self.tmp_annot_filters]

        # Execute command
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = list()
        for record in self.variants:
            if record.info["is_filtered"] == 0:
                expected.append(record.id)
                if "ANN" in record.info:
                    for curr_ann in record.info["ANN"]:
                        if curr_ann["is_filtered"] == 0:
                            expected.append(curr_ann["id"])
        observed = list()
        with AnnotVCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
                if "ANN" in record.info:
                    for curr_ann in record.info["ANN"]:
                        observed.append(curr_ann["id"])
        self.assertEqual(observed, expected)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
