#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import json
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
class TestCreateModel(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_out_info = os.path.join(tmp_folder, unique_id + "_out_info.tsv")
        self.tmp_out_model = os.path.join(tmp_folder, unique_id + "_out_model.json")
        # Lengths distributions
        self.tmp_len = os.path.join(tmp_folder, unique_id + "_len.json")
        with open(self.tmp_len, "w") as writer:
            writer.write("""[
{
    "loci": {
        "11:102322740-102322841": {
            "name": "NR27",
            "position": "11:102322740-102322841",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"91": 2, "92": 2, "93": 2, "94": 3, "95": 1, "96": 3, "99": 1},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        },
        "2:47414383-47414484": {
            "name": "BAT26",
            "position": "2:47414383-47414484",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"91": 2, "92": 2, "93": 2, "94": 4, "95": 6, "96": 1, "97": 2, "98": 3},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        },
        "4:54732007-54732108": {
            "name": "BAT25",
            "position": "4:54732007-54732108",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"92": 2, "94": 1, "95": 1, "96": 2, "97": 1, "98": 3, "99": 12, "100": 9, "101": 1, "102": 3},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        }
    },
    "name": "splA",
    "results": {}
},
{
    "loci": {
        "11:102322740-102322841": {
            "name": "NR27",
            "position": "11:102322740-102322841",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"83": 2, "84": 4, "85": 21, "86": 36, "87": 21, "88": 34, "89": 44, "90": 4, "93": 1, "94": 3, "95": 3, "96": 3},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        },
        "2:47414383-47414484": {
            "name": "BAT26",
            "position": "2:47414383-47414484",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"82": 7, "83": 14, "84": 49, "85": 100, "86": 61, "87": 64, "88": 81, "89": 14, "91": 1, "92": 2, "95": 1, "98": 1},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        },
        "4:54732007-54732108": {
            "name": "BAT25",
            "position": "4:54732007-54732108",
            "results": {
                "model": {
                    "data": {
                        "lengths": {
                            "ct_by_len": {"88": 3, "89": 2, "90": 5, "91": 21, "92": 15, "93": 9, "94": 8, "95": 6, "96": 3, "97": 11, "98": 3, "99": 8, "100": 2, "101": 4, "102": 2, "103": 3},
                            "mode": "reads"
                        }
                    },
                    "score": null,
                    "status": null
                }
            }
        }
    },
    "name": "splB",
    "results": {}
}
]""")
        # Microsatellites
        self.tmp_ms = os.path.join(tmp_folder, unique_id + "_ms.bed")
        with open(self.tmp_ms, "w") as writer:
            writer.write("""2	47414383	47414484	BAT26
4	54732007	54732108	BAT25
11	102322740	102322841	NR27""")
        # Status
        self.tmp_status = os.path.join(tmp_folder, unique_id + "_status.tsv")
        with open(self.tmp_status, "w") as writer:
            writer.write("""sample	locus_position	method_id	key	value	type
splA	2:47414383-47414484	model	status	MSS	str
splA	4:54732007-54732108	model	status	MSS	str
splA	11:102322740-102322841	model	status	MSS	str
splB	2:47414383-47414484	model	status	MSI	str
splB	4:54732007-54732108	model	status	MSI	str
splB	11:102322740-102322841	model	status	MSI	str""")

    def tearDown(self):
        for curr_file in [self.tmp_out_info, self.tmp_out_model, self.tmp_len, self.tmp_ms, self.tmp_status]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_all(self):
        cmd = [
            "microsatCreateModel.py",
            "--min-support", str(10),
            "--peak-height-cutoff", str(0.05),
            "--inputs-length-distributions", self.tmp_len,
            "--input-loci-status", self.tmp_status,
            "--input-microsatellites", self.tmp_ms,
            "--output-info", self.tmp_out_info,
            "--output-model", self.tmp_out_model
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        # Info
        expected_info = """Nb retained samples: 2
Locus_position	Locus_name	MSI	MSS	Unused
11:102322740-102322841	NR27	1	1	0
2:47414383-47414484	BAT26	1	1	0
4:54732007-54732108	BAT25	1	1	0""".split("\n")
        with open(self.tmp_out_info) as reader:
            observed_info = [elt.strip() for elt in reader.readlines()]
        self.assertEqual(observed_info, expected_info)
        # Model
        expected_model = [
            {
                "loci": {
                    "11:102322740-102322841": {
                        "name": "NR27",
                        "position": "11:102322740-102322841",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.069307, "pro_q": 0.0},
                                    "lengths": {
                                        "ct_by_len": {91: 2, 92: 2, 93: 2, 94: 3, 95: 1, 96: 3, 99: 1},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 7, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSS"
                            }
                        }
                    },
                    "2:47414383-47414484": {
                        "name": "BAT26",
                        "position": "2:47414383-47414484",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.063006, "pro_q": 0.0},
                                    "lengths": {
                                        "ct_by_len": {91: 2, 92: 2, 93: 2, 94: 4, 95: 6, 96: 1, 97: 2, 98: 3},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 8, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSS"
                            }
                        }
                    },
                    "4:54732007-54732108": {
                        "name": "BAT25",
                        "position": "4:54732007-54732108",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.02459, "pro_q": 0.000848},
                                    "lengths": {
                                        "ct_by_len": {92: 2, 94: 1, 95: 1, 96: 2, 97: 1, 98: 3, 99: 12, 100: 9, 101: 1, 102: 3},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSS"
                            }
                        }
                    }
                },
                "name": "splA",
                "results": {}
            },
            {
                "loci": {
                    "11:102322740-102322841": {
                        "name": "NR27",
                        "position": "11:102322740-102322841",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.132201, "pro_q": 0.0},
                                    "lengths": {
                                        "ct_by_len": {83: 2, 84: 4, 85: 21, 86: 36, 87: 21, 88: 34, 89: 44, 90: 4, 93: 1, 94: 3, 95: 3, 96: 3},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSI"
                            }
                        }
                    },
                    "2:47414383-47414484": {
                        "name": "BAT26",
                        "position": "2:47414383-47414484",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.147562, "pro_q": 0.0},
                                    "lengths": {
                                        "ct_by_len": {82: 7, 83: 14, 84: 49, 85: 100, 86: 61, 87: 64, 88: 81, 89: 14, 91: 1, 92: 2, 95: 1, 98: 1},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 8, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSI"
                            }
                        }
                    },
                    "4:54732007-54732108": {
                        "name": "BAT25",
                        "position": "4:54732007-54732108",
                        "results": {
                            "model": {
                                "data": {
                                    "MSIsensor-pro": {"pro_p": 0.067182, "pro_q": 0.000754},
                                    "lengths": {
                                        "ct_by_len": {88: 3, 89: 2, 90: 5, 91: 21, 92: 15, 93: 9, 94: 8, 95: 6, 96: 3, 97: 11, 98: 3, 99: 8, 100: 2, 101: 4, 102: 2, 103: 3},
                                        "mode": "reads"
                                    },
                                    "mSINGS": {"nb_peaks": 16, "peak_height_cutoff": 0.05}
                                },
                                "score": None,
                                "status": "MSI"
                            }
                        }
                    }
                },
                "name": "splB",
                "results": {}
            }
        ]
        with open(self.tmp_out_model) as reader:
            observed_model = json.load(reader)
            for spl in observed_model:
                for locus_id, locus in spl["loci"].items():
                    data = locus["results"]["model"]["data"]
                    data["MSIsensor-pro"]["pro_p"] = round(data["MSIsensor-pro"]["pro_p"], 6)
                    data["MSIsensor-pro"]["pro_q"] = round(data["MSIsensor-pro"]["pro_q"], 6)
                    data["lengths"]["ct_by_len"] = {int(len): ct for len, ct in data["lengths"]["ct_by_len"].items()}
        for curr_obs, curr_expect in zip(observed_model, expected_model):
            self.assertDictEqual(curr_obs, curr_expect)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
