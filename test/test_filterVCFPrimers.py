#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess
from anacore.bed import BEDIO, BEDRecord
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr
from anacore.sequenceIO import FastaIO, Sequence

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class FilterVCFPrimers(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_sequences = os.path.join(tmp_folder, unique_id + ".fasta")
        self.tmp_regions = os.path.join(tmp_folder, unique_id + ".bed")
        self.tmp_variants = os.path.join(tmp_folder, unique_id + ".vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "filterVCFPrimers.py",
            "--input-variants", self.tmp_variants,
            "--input-regions", self.tmp_regions,
            "--input-sequences", self.tmp_sequences,
            "--output-variants", self.tmp_output
        ]

        # Create fasta
        with FastaIO(self.tmp_sequences, "w") as FH_seq:
            FH_seq.write(Sequence("artificial_chr1", "NNNAAAATTTGGGGGGGGGGTTTAAANNN"))
            #                                         123456789| | | | | | | | | |
            #                                                  10| 14| 18| 22| 26|
            #                                                    12  16  20  24  28
            FH_seq.write(Sequence("artificial_chr2", "CGATNNNCGAT"))
            #                                         123456789|
            #                                                  10

        # Create VCF
        with VCFIO(self.tmp_variants, "w") as FH_var:
            FH_var.info = {"ZOI": HeaderInfoAttr("ZOI", "If the variant can be in interest area.", type="String", number="1")}
            FH_var.writeHeader()
            self.variants = [
                VCFRecord("artificial_chr1", 6, "alt_0", "A", ["AA"], None, None, {"ZOI": "no"}),
                VCFRecord("artificial_chr1", 8, "alt_1", "TT", ["T"], None, None, {"ZOI": "no"}),
                VCFRecord("artificial_chr1", 8, "alt_2", "T", ["TT"], None, None, {"ZOI": "yes"}),
                VCFRecord("artificial_chr1", 9, "alt_3", "TTGG", ["TT"], None, None, {"ZOI": "yes"}),
                VCFRecord("artificial_chr1", 14, "alt_4", "G", ["GG"], None, None, {"ZOI": "yes"}),
                VCFRecord("artificial_chr1", 18, "alt_5", "GGG", ["G"], None, None, {"ZOI": "yes"}),  # ZOI downstream limit deletion
                VCFRecord("artificial_chr1", 22, "alt_6", "T", ["TT"], None, None, {"ZOI": "yes"}),

                VCFRecord("artificial_chr1", 9, "alt_7", "TT", ["TC"], None, None, {"ZOI": "no"}),  # Substitution before end of upstream primer
                VCFRecord("artificial_chr1", 10, "alt_8", "TG", ["TC"], None, None, {"ZOI": "yes"}),  # Substitution in upstream limit of ZOI
                VCFRecord("artificial_chr1", 15, "alt_9", "GG", ["GC"], None, None, {"ZOI": "yes"}),  # Substitution in dosnstream limit of ZOI
                VCFRecord("artificial_chr1", 20, "alt_10", "GT", ["GC"], None, None, {"ZOI": "no"}),  # Substitution after start of downstream primer
                VCFRecord("artificial_chr1", 21, "alt_11", "TT", ["TC"], None, None, {"ZOI": "no"}),  # Substitution in downstream primer

                VCFRecord("artificial_chr2", 1, "alt_12", "C", ["CTT"], None, None, {"ZOI": "no"}),  # Insertion before end of upstream primer
                VCFRecord("artificial_chr2", 2, "alt_13", "G", ["GCC"], None, None, {"ZOI": "yes"}),  # Insertion in upstream limit of ZOI
                VCFRecord("artificial_chr2", 3, "alt_14", "AT", ["CCGC"], None, None, {"ZOI": "yes"}),  # Insertion in upstream limit of ZOI and without standardization
                VCFRecord("artificial_chr2", 9, "alt_15", "G", ["GCC"], None, None, {"ZOI": "yes"}),  # Insertion in downstream limit of ZOI
                VCFRecord("artificial_chr2", 9, "alt_16", "G", ["NNN"], None, None, {"ZOI": "yes"}),  # Insertion in downstream limit of ZOI and without standardization
                VCFRecord("artificial_chr2", 10, "alt_17", "-", ["CC"], None, None, {"ZOI": "yes"}),  # Insertion in downstream limit of ZOI
                VCFRecord("artificial_chr2", 10, "alt_18", "A", ["ATT"], None, None, {"ZOI": "no"}),  # Insertion after start of downstream primer

                VCFRecord("artificial_chr2", 1, "alt_19", "CG", ["C"], None, None, {"ZOI": "no"}),  # Deletion before end of upstream primer
                VCFRecord("artificial_chr2", 2, "alt_20", "GA", ["G"], None, None, {"ZOI": "yes"}),  # Deletion in upstream limit of ZOI
                VCFRecord("artificial_chr2", 3, "alt_21", "AT", ["C"], None, None, {"ZOI": "yes"}),  # Deletion in upstream limit of ZOI and without standardization
                VCFRecord("artificial_chr2", 6, "alt_22", "NNCG", ["N"], None, None, {"ZOI": "yes"}),  # Deletion in downstream limit of ZOI
                VCFRecord("artificial_chr2", 8, "alt_23", "CG", ["C"], None, None, {"ZOI": "yes"}),  # Deletion in downstream limit of ZOI
                VCFRecord("artificial_chr2", 8, "alt_24", "CG", ["T"], None, None, {"ZOI": "yes"}),  # Deletion in downstream limit of ZOI and without standardization
                VCFRecord("artificial_chr2", 9, "alt_25", "GA", ["G"], None, None, {"ZOI": "no"}),  # Insertion after start of downstream primer
                VCFRecord("artificial_chr2", 10, "alt_26", "A", ["-"], None, None, {"ZOI": "no"}),  # Insertion after start of downstream primer
                VCFRecord("artificial_chr2", 10, "alt_27", "AT", ["A"], None, None, {"ZOI": "no"}),  # Insertion after start of downstream primer
            ]
            for idx, curr_var in enumerate(self.variants):
                FH_var.write(curr_var)

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_sequences, self.tmp_regions, self.tmp_variants, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def testOverlapException(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord("artificial_chr1", 5, 25, "ampl1", None, "+", 11, 20)
            FH_reg.write(ampl1)
            ampl2 = BEDRecord("artificial_chr2", 1, 11, "ampl2", None, "+", 3, 9)
            FH_reg.write(ampl2)
            ampl3 = BEDRecord("artificial_chr1", 23, 28, "ampl2", None, "+", 25, 26)
            FH_reg.write(ampl3)

        # Execute command
        with self.assertRaises(subprocess.CalledProcessError) as context:
            subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

    def testResults(self):
        # Create BED
        with BEDIO(self.tmp_regions, "w", 8) as FH_reg:
            ampl1 = BEDRecord("artificial_chr1", 5, 25, "ampl1", None, "+", 11, 20)
            FH_reg.write(ampl1)
            ampl2 = BEDRecord("artificial_chr2", 1, 11, "ampl2", None, "+", 3, 9)
            FH_reg.write(ampl2)

        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = [curr_var.id for curr_var in self.variants if curr_var.info["ZOI"] == "yes"]
        observed = list()
        with VCFIO(self.tmp_output) as FH_results:
            for record in FH_results:
                observed.append(record.id)
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
