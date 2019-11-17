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
from anacore.vcf import HeaderInfoAttr, VCFIO, VCFRecord

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class TestUniformCallersHaplotypes(unittest.TestCase):
    def tearDown(self):
        # Clean temporary files
        for curr_caller in sorted(self.callers):
            for curr_pattern in [self.tmp_initial_pathes, self.tmp_haplotyped_pathes, self.tmp_expected_pathes, self.tmp_out_pathes]:
                curr_file = curr_pattern.format(curr_caller)
                if os.path.exists(curr_file):
                    os.remove(curr_file)

    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_initial_pathes = os.path.join(tmp_folder, unique_id + "_{}_initial.vcf")
        self.tmp_haplotyped_pathes = os.path.join(tmp_folder, unique_id + "_{}_haplotyped.vcf")
        self.tmp_expected_pathes = os.path.join(tmp_folder, unique_id + "_{}_expected.vcf")
        self.tmp_out_pathes = os.path.join(tmp_folder, unique_id + "_{}_out.vcf")

        # test cases
        self.test_cases = [
            {  # *a-b, a-b, a b, /
                "initial": {
                    "caller1": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller2": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller3": [
                        VCFRecord("chr1", 14, None, "G", ["C"], info={"AD": 100}),
                        VCFRecord("chr1", 18, None, "A", ["G"], info={"AD": 104})
                    ]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller2": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller3": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"], info={"MCO_VAR": ["chr1:14=G/C", "chr1:18=A/G"], "AD": 100})]
                },
                "expected": {
                    "caller1": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller2": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"])],
                    "caller3": [VCFRecord("chr1", 14, None, "GCGTA", ["CCGTG"], info={"AD": 104})]
                }
            },
            {  # *a b, a b, a-b, /
                "initial": {
                    "caller1": [
                        VCFRecord("chr2", 14, None, "G", ["C"]),
                        VCFRecord("chr2", 18, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr2", 14, None, "G", ["C"]),
                        VCFRecord("chr2", 18, None, "A", ["G"])
                    ],
                    "caller3": [VCFRecord("chr2", 14, None, "GCGTA", ["CCGTG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr2", 14, None, "GCGTA", ["CCGTG"], info={"MCO_VAR": ["chr2:14=G/C", "chr2:18=A/G"]})],
                    "caller2": [VCFRecord("chr2", 14, None, "GCGTA", ["CCGTG"], info={"MCO_VAR": ["chr2:14=G/C", "chr2:18=A/G"]})],
                    "caller3": [VCFRecord("chr2", 14, None, "GCGTA", ["CCGTG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr2", 14, None, "G", ["C"]),
                        VCFRecord("chr2", 18, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr2", 14, None, "G", ["C"]),
                        VCFRecord("chr2", 18, None, "A", ["G"])
                    ],
                    "caller3": [
                        VCFRecord("chr2", 14, None, "G", ["C"]),
                        VCFRecord("chr2", 18, None, "A", ["G"])
                    ]
                }
            },
            {  # *a-b c, a-b c, a b c, /
                "initial": {
                    "caller1": [
                        VCFRecord("chr3", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr3", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr3", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr3", 20, None, "A", ["G"])
                    ],
                    "caller3": [
                        VCFRecord("chr3", 14, None, "G", ["C"], info={"AD": 104}),
                        VCFRecord("chr3", 18, None, "A", ["G"], info={"AD": 100}),
                        VCFRecord("chr3", 20, None, "A", ["G"], info={"AD": 98})
                    ]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr3", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr3:14=GCGTA/CCGTG", "chr3:20=A/G"]})],
                    "caller2": [VCFRecord("chr3", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr3:14=GCGTA/CCGTG", "chr3:20=A/G"]})],
                    "caller3": [VCFRecord("chr3", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr3:14=G/C", "chr3:18=A/G", "chr3:20=A/G"], "AD": 98})]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr3", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr3", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr3", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr3", 20, None, "A", ["G"])
                    ],
                    "caller3": [
                        VCFRecord("chr3", 14, None, "GCGTA", ["CCGTG"], info={"AD": 104}),
                        VCFRecord("chr3", 20, None, "A", ["G"], info={"AD": 98})
                    ]
                }
            },
            {  # *a-b c, a-b c, a b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr4", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr4", 20, None, "A", ["G"])
                    ],
                    "caller3": [
                        VCFRecord("chr4", 14, None, "G", ["C"], info={"AD": 98}),
                        VCFRecord("chr4", 18, None, "A", ["G"], info={"AD": 104}),
                        VCFRecord("chr4", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [VCFRecord("chr4", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr4", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr4:14=GCGTA/CCGTG", "chr4:20=A/G"]})],
                    "caller2": [VCFRecord("chr4", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr4:14=GCGTA/CCGTG", "chr4:20=A/G"]})],
                    "caller3": [VCFRecord("chr4", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr4:14=G/C", "chr4:18=A/G", "chr4:20=A/G"], "AD": 98})],
                    "caller4": [VCFRecord("chr4", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr4", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr4", 20, None, "A", ["G"])
                    ],
                    "caller3": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"], info={"AD": 104}),
                        VCFRecord("chr4", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [
                        VCFRecord("chr4", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr4", 20, None, "A", ["G"])
                    ]
                }
            },
            {  # *a-b c, a' a-b c, a b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr5", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr5", 14, None, "G", ["C"], info={"AD": 3}),
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr5", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr5", 14, None, "G", ["C"], info={"AD": 110}),
                        VCFRecord("chr5", 18, None, "A", ["G"], info={"AD": 105}),
                        VCFRecord("chr5", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [VCFRecord("chr5", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr5", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr5:14=GCGTA/CCGTG", "chr5:20=A/G"]})],
                    "caller2": [
                        VCFRecord("chr5", 14, None, "G", ["C"], info={"AD": 3}),
                        VCFRecord("chr5", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr5:14=GCGTA/CCGTG", "chr5:20=A/G"], "AD": 100})
                    ],
                    "caller3": [VCFRecord("chr5", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr5:14=G/C", "chr5:18=A/G", "chr5:20=A/G"], "AD": 100})],
                    "caller4": [VCFRecord("chr5", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr5", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr5", 14, None, "G", ["C"], info={"AD": 3}),
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr5", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"], info={"AD": 110}),
                        VCFRecord("chr5", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [
                        VCFRecord("chr5", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr5", 20, None, "A", ["G"])
                    ]
                }
            },
            {  # *a b c, a' a-b c, a-b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr6", 14, None, "G", ["C"]),
                        VCFRecord("chr6", 18, None, "A", ["G"]),
                        VCFRecord("chr6", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr6", 14, None, "G", ["C"], info={"AD": 3}),
                        VCFRecord("chr6", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr6", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr6", 14, None, "GCGTA", ["CCGTG"], info={"AD": 105}),
                        VCFRecord("chr6", 20, None, "A", ["G"], info={"AD": 101})
                    ],
                    "caller4": [VCFRecord("chr6", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr6", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr6:14=G/C", "chr6:18=A/G", "chr6:20=A/G"]})],
                    "caller2": [
                        VCFRecord("chr6", 14, None, "G", ["C"], info={"AD": 3}),
                        VCFRecord("chr6", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr6:14=GCGTA/CCGTG", "chr6:20=A/G"], "AD": 100})
                    ],
                    "caller3": [VCFRecord("chr6", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr6:14=GCGTA/CCGTG", "chr6:20=A/G"], "AD": 101})],
                    "caller4": [VCFRecord("chr6", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr6", 14, None, "G", ["C"]),
                        VCFRecord("chr6", 18, None, "A", ["G"]),
                        VCFRecord("chr6", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr6", 14, None, "G", ["C"], info={"AD": 100}),
                        VCFRecord("chr6", 18, None, "A", ["G"], info={"AD": 100}),
                        VCFRecord("chr6", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr6", 14, None, "G", ["C"], info={"AD": 105}),
                        VCFRecord("chr6", 18, None, "A", ["G"], info={"AD": 105}),
                        VCFRecord("chr6", 20, None, "A", ["G"], info={"AD": 101})
                    ],
                    "caller4": [
                        VCFRecord("chr6", 14, None, "G", ["C"]),
                        VCFRecord("chr6", 18, None, "A", ["G"]),
                        VCFRecord("chr6", 20, None, "A", ["G"])
                    ]
                }
            },
            {  # *a b c, a-b b' c, a-b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr7", 14, None, "G", ["C"]),
                        VCFRecord("chr7", 18, None, "A", ["G"]),
                        VCFRecord("chr7", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr7", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr7", 18, None, "A", ["G"], info={"AD": 3}),
                        VCFRecord("chr7", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr7", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr7", 20, None, "A", ["G"])
                    ],
                    "caller4": [VCFRecord("chr7", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr7", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr7:14=G/C", "chr7:18=A/G", "chr7:20=A/G"]})],
                    "caller2": [
                        VCFRecord("chr7", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr7:14=GCGTA/CCGTG", "chr7:20=A/G"], "AD": 100}),
                        VCFRecord("chr7", 18, None, "G", ["C"], info={"AD": 3})
                    ],
                    "caller3": [VCFRecord("chr7", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr7:14=GCGTA/CCGTG", "chr7:20=A/G"]})],
                    "caller4": [VCFRecord("chr7", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr7", 14, None, "G", ["C"]),
                        VCFRecord("chr7", 18, None, "A", ["G"]),
                        VCFRecord("chr7", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr7", 14, None, "G", ["C"], info={"AD": 100}),
                        VCFRecord("chr7", 18, None, "A", ["G"], info={"AD": 100}),
                        VCFRecord("chr7", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr7", 14, None, "G", ["C"]),
                        VCFRecord("chr7", 18, None, "A", ["G"]),
                        VCFRecord("chr7", 20, None, "A", ["G"])
                    ],
                    "caller4": [
                        VCFRecord("chr7", 14, None, "G", ["C"]),
                        VCFRecord("chr7", 18, None, "A", ["G"]),
                        VCFRecord("chr7", 20, None, "A", ["G"])
                    ]
                }
            },
            {  # *a-b c, a-b b' c, a b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr8", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr8", 18, None, "A", ["G"], info={"AD": 3}),
                        VCFRecord("chr8", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr8", 14, None, "G", ["C"], info={"AD": 110}),
                        VCFRecord("chr8", 18, None, "A", ["G"], info={"AD": 105}),
                        VCFRecord("chr8", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [VCFRecord("chr8", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [VCFRecord("chr8", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr8:14=GCGTA/CCGTG", "chr8:20=A/G"]})],
                    "caller2": [
                        VCFRecord("chr8", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr8:14=GCGTA/CCGTG", "chr8:20=A/G"], "AD": 100}),
                        VCFRecord("chr8", 18, None, "G", ["C"], info={"AD": 3})
                    ],
                    "caller3": [VCFRecord("chr8", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr8:14=G/C", "chr8:18=A/G", "chr8:20=A/G"], "AD": 100})],
                    "caller4": [VCFRecord("chr8", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr8", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr8", 18, None, "A", ["G"], info={"AD": 3}),
                        VCFRecord("chr8", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"], info={"AD": 110}),
                        VCFRecord("chr8", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [
                        VCFRecord("chr8", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr8", 20, None, "A", ["G"])
                    ]
                }
            },
            {  # *a' a-b c, a-b b' c, a b c, a-b-c
                "initial": {
                    "caller1": [
                        VCFRecord("chr9", 14, None, "G", ["C"]),
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr9", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr9", 18, None, "A", ["G"], info={"AD": 3}),
                        VCFRecord("chr9", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr9", 14, None, "G", ["C"], info={"AD": 110}),
                        VCFRecord("chr9", 18, None, "A", ["G"], info={"AD": 105}),
                        VCFRecord("chr9", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [VCFRecord("chr9", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "haplotyped": {
                    "caller1": [
                        VCFRecord("chr9", 14, None, "G", ["C"]),
                        VCFRecord("chr9", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr9:14=GCGTA/CCGTG", "chr9:20=A/G"]})
                    ],
                    "caller2": [
                        VCFRecord("chr9", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr9:14=GCGTA/CCGTG", "chr9:20=A/G"], "AD": 100}),
                        VCFRecord("chr9", 18, None, "G", ["C"], info={"AD": 3})
                    ],
                    "caller3": [VCFRecord("chr9", 14, None, "GCGTATCA", ["CCGTGTCG"], info={"MCO_VAR": ["chr9:14=G/C", "chr9:18=A/G", "chr9:20=A/G"], "AD": 100})],
                    "caller4": [VCFRecord("chr9", 14, None, "GCGTATCA", ["CCGTGTCG"])]
                },
                "expected": {
                    "caller1": [
                        VCFRecord("chr9", 14, None, "G", ["C"]),
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr9", 20, None, "A", ["G"])
                    ],
                    "caller2": [
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"], info={"AD": 100}),
                        VCFRecord("chr9", 18, None, "A", ["G"], info={"AD": 3}),
                        VCFRecord("chr9", 20, None, "A", ["G"], info={"AD": 104})
                    ],
                    "caller3": [
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"], info={"AD": 110}),
                        VCFRecord("chr9", 20, None, "A", ["G"], info={"AD": 100})
                    ],
                    "caller4": [
                        VCFRecord("chr9", 14, None, "GCGTA", ["CCGTG"]),
                        VCFRecord("chr9", 20, None, "A", ["G"])
                    ]
                }
            }
        ]

        # Get callers
        callers = set()
        for curr_test in self.test_cases:
            for curr_caller in curr_test["initial"]:
                callers.add(curr_caller)
        self.callers = sorted(list(callers))

        # Write files
        for curr_caller in self.callers:
            # Initial
            with VCFIO(self.tmp_initial_pathes.format(curr_caller), "w") as handle_out:
                handle_out.info = {
                    "AD": HeaderInfoAttr("AD", "Alternative allele depth.", type="Integer", number="1")
                }
                handle_out.extra_header = ["##source={}".format(curr_caller)]
                handle_out.writeHeader()
                for curr_test in self.test_cases:
                    if curr_caller in curr_test["initial"]:
                        for curr_var in curr_test["initial"][curr_caller]:
                            handle_out.write(curr_var)
            # Haplotyped
            with VCFIO(self.tmp_haplotyped_pathes.format(curr_caller), "w") as handle_out:
                handle_out.info = {
                    "AD": HeaderInfoAttr("AD", "Alternative allele depth.", type="Integer", number="1"),
                    "MCO_VAR": HeaderInfoAttr("MCO_VAR", "Name of the variants merged because their occur on same reads.", type="String", number=".")
                }
                handle_out.extra_header = ["##source={}".format(curr_caller)]
                handle_out.writeHeader()
                for curr_test in self.test_cases:
                    if curr_caller in curr_test["haplotyped"]:
                        for curr_var in curr_test["haplotyped"][curr_caller]:
                            handle_out.write(curr_var)
            # Expected
            with VCFIO(self.tmp_expected_pathes.format(curr_caller), "w") as handle_out:
                handle_out.info = {
                    "AD": HeaderInfoAttr("AD", "Alternative allele depth.", type="Integer", number="1")
                }
                handle_out.extra_header = ["##source={}".format(curr_caller)]
                handle_out.writeHeader()
                for curr_test in self.test_cases:
                    if curr_caller in curr_test["expected"]:
                        for curr_var in curr_test["expected"][curr_caller]:
                            handle_out.write(curr_var)

    def testUniformCallersHaplotypes(self):
        # Execute command
        cmd = [
            "uniformCallersHaplotypes.py",
            "--calling-sources", *self.callers,
            "--inputs-haplotyped", *[self.tmp_haplotyped_pathes.format(curr_caller) for curr_caller in self.callers],
            "--inputs-variants", *[self.tmp_initial_pathes.format(curr_caller) for curr_caller in self.callers],
            "--outputs-variants", *[self.tmp_out_pathes.format(curr_caller) for curr_caller in self.callers]
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)

        # Validate results
        for curr_caller in self.callers:
            expected = []
            with open(self.tmp_expected_pathes.format(curr_caller)) as handle_in:
                expected = [elt.replace("\t\t", "\t.\t") for elt in handle_in.readlines()]
            observed = []
            with open(self.tmp_out_pathes.format(curr_caller)) as handle_in:
                observed = handle_in.readlines()
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
