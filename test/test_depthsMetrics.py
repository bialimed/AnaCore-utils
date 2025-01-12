#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import os
import sys
import tempfile
import unittest
import uuid

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from depthsMetrics import getDistribution, loadFromDepthFile


########################################################################
#
# FUNCTIONS
#
########################################################################
class DepthsMetrics(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_dp = os.path.join(tmp_folder, unique_id + "_depths.tsv")
        with open(self.tmp_dp, "w") as writer:
            writer.write("""1	100	0	10
1	101	0	10
1	102	0	10
1	103	5	50
1	104	6	60
1	105	4	40
1	106	7	70
1	107	5	50
1	108	5	50
1	109	6	60
1	110	4	40
2	210	12	120
2	211	13	13
2	212	11	110
2	213	12	120
2	214	11	110
2	215	12	120
2	216	10	100
2	217	11	110
2	218	12	0
2	219	13	0
2	220	15	150""")

    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_dp]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_getDistribution(self):
        datasets = [
            {
                "data": {x: 1 for x in range(1, 11)},
                "expected": {
                    10: {'min': 1, 'max': 10, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 3.5, '40_percentile': 4.5, '50_percentile': 5.5, '60_percentile': 6.5, '70_percentile': 7.5, '80_percentile': 8.5, '90_percentile': 9.5},
                    25: {'min': 1, 'max': 10, '25_percentile': 3, '50_percentile': 5.5, '75_percentile': 8},
                    33.3: {'min': 1, 'max': 10, '33.3_percentile': 3.5, '66.6_percentile': 7.5},
                }
            },
            {
                "data": {x: 2 for x in range(1, 11)},
                "expected": {
                    10: {'min': 1, 'max': 10, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 3.5, '40_percentile': 4.5, '50_percentile': 5.5, '60_percentile': 6.5, '70_percentile': 7.5, '80_percentile': 8.5, '90_percentile': 9.5},
                    25: {'min': 1, 'max': 10, '25_percentile': 3, '50_percentile': 5.5, '75_percentile': 8},
                    33.3: {'min': 1, 'max': 10, '33.3_percentile': 4, '66.6_percentile': 7},
                }
            },
            {
                "data": {x: 3 for x in range(1, 11)},
                "expected": {
                    10: {'min': 1, 'max': 10, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 3.5, '40_percentile': 4.5, '50_percentile': 5.5, '60_percentile': 6.5, '70_percentile': 7.5, '80_percentile': 8.5, '90_percentile': 9.5},
                    25: {'min': 1, 'max': 10, '25_percentile': 3, '50_percentile': 5.5, '75_percentile': 8},
                    33.3: {'min': 1, 'max': 10, '33.3_percentile': 4, '66.6_percentile': 7},
                }
            },
            {
                "data": {x: 4 for x in range(1, 11)},
                "expected": {
                    10: {'min': 1, 'max': 10, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 3.5, '40_percentile': 4.5, '50_percentile': 5.5, '60_percentile': 6.5, '70_percentile': 7.5, '80_percentile': 8.5, '90_percentile': 9.5},
                    25: {'min': 1, 'max': 10, '25_percentile': 3, '50_percentile': 5.5, '75_percentile': 8},
                    33.3: {'min': 1, 'max': 10, '33.3_percentile': 4, '66.6_percentile': 7},
                }
            },
            {
                "data": {x: 1 for x in [11, 1, 3, 6, 2, 10, 7, 8, 9]},
                "expected": {
                    10: {'min': 1, 'max': 11, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 4.5, '40_percentile': 6.5, '50_percentile': 7, '60_percentile': 7.5, '70_percentile': 8.5, '80_percentile': 9.5, '90_percentile': 10.5},
                    25: {'min': 1, 'max': 11, '25_percentile': 2.5, '50_percentile': 7, '75_percentile': 9.5},
                    33.3: {'min': 1, 'max': 11, '33.3_percentile': 4.5, '66.6_percentile': 8.5},
                }
            },
            {
                "data": {1: 1, 2: 1, 3: 2, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1},
                "expected": {
                    10: {'min': 1, 'max': 11, '10_percentile': 1.5, '20_percentile': 2.5, '30_percentile': 3, '40_percentile': 4.5, '50_percentile': 6.5, '60_percentile': 7.5, '70_percentile': 8.5, '80_percentile': 9.5, '90_percentile': 10.5},
                    25: {'min': 1, 'max': 11, '25_percentile': 3, '50_percentile': 6.5, '75_percentile': 9},
                    33.3: {'min': 1, 'max': 11, '33.3_percentile': 3, '66.6_percentile': 8.5},
                }
            }
        ]
        for curr in datasets:
            for percentile_step in [10, 25, 33.3]:
                self.assertEqual(
                    getDistribution(curr["data"], percentile_step, 4),
                    curr["expected"][percentile_step]
                )

    def test_loadFromDepthFile(self):
        dp_list, ct_list_by_spl = loadFromDepthFile(self.tmp_dp, ["splA", "splB"])
        self.assertEqual(
            dp_list,
            [0, 4, 5, 6, 7, 10, 11, 12, 13, 15, 40, 50, 60, 70, 100, 110, 120, 150]
        )
        self.assertEqual(
            ct_list_by_spl,
            {
                "splA": [3, 2, 3, 2, 1, 1, 3, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                "splB": [2, 0, 0, 0, 0, 3, 0, 0, 1, 0, 2, 3, 2, 1, 1, 3, 3, 1]
            }
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
