#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2021 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'

import numpy
import os
import statistics
import sys
import unittest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
APP_DIR = os.path.dirname(TEST_DIR)
BIN_DIR = os.path.join(APP_DIR, "bin")
sys.path.append(BIN_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

from depthsMetrics import getDistribution


########################################################################
#
# FUNCTIONS
#
########################################################################
def getDistributionNumpy(ct_by_dp, percentile_step=25, precision=4):
    values = []
    for dp, count in ct_by_dp.items():
        for idx in range(count):
            values.append(dp)
    distrib = {
        "min": round(min(values), precision),
        "max": round(max(values), precision)
    }
    for curr_percentile in range(percentile_step, 100, percentile_step):
        distrib['{:02}'.format(curr_percentile) + "_percentile"] = round(numpy.percentile(values, curr_percentile, interpolation="midpoint"), precision)
    return distrib


def getDistributionStatistics(ct_by_dp, percentile_step=25, precision=4):
    values = []
    for dp, count in ct_by_dp.items():
        for idx in range(count):
            values.append(dp)
    distrib = {
        "min": round(min(values), precision),
        "max": round(max(values), precision)
    }
    values = statistics.quantiles(values, n=100, method='inclusive')
    for curr_percentile in range(percentile_step, 100, percentile_step):
        distrib['{:02}'.format(curr_percentile) + "_percentile"] = round(values[curr_percentile + 1], precision)
    return distrib


class DepthsMetrics(unittest.TestCase):
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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
