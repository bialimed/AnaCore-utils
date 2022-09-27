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
        self.tmp_in_msings = os.path.join(tmp_folder, unique_id + "_in_msings.json")
        self.tmp_in_msisensor = os.path.join(tmp_folder, unique_id + "_in_msisensor.json")
        self.tmp_in_random = os.path.join(tmp_folder, unique_id + "_in_random.json")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.json")
        with open(self.tmp_in_msings, "w") as writer:
            writer.write("""[
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "mSINGS": {"data": {}, "score": null, "status": "MSS"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "mSINGS": {"data": {}, "score": null, "status": "MSS"}
                }
            }
        },
        "name": "splA",
        "results": {
            "mSINGS": {"score": 1.0, "status": "MSS"}
        }
    },
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "mSINGS": {"data": {}, "score": null, "status": "MSS"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "mSINGS": {"data": {}, "score": null, "status": "MSI"}
                }
            }
        },
        "name": "splB",
        "results": {
            "mSINGS": {"score": 0.5, "status": "MSI"}
        }
    }
]""")
        with open(self.tmp_in_msisensor, "w") as writer:
            writer.write("""[
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "MSISensor": {"data": {}, "score": 1.0, "status": "MSS"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "MSISensor": {"data": {}, "score": 1.0, "status": "MSS"}
                }
            }
        },
        "name": "splA",
        "results": {
            "MSISensor": {"score": 1.0, "status": "MSS"}
        }
    },
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "MSISensor": {"data": {}, "score": 0.9, "status": "MSI"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "MSISensor": {"data": {}, "score": 0.8, "status": "MSI"}
                }
            }
        },
        "name": "splB",
        "results": {
            "MSISensor": {"score": 0.95, "status": "MSI"}
        }
    }
]""")
        with open(self.tmp_in_random, "w") as writer:
            writer.write("""[
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "RandomForest": {"data": {}, "score": 0.99, "status": "MSS"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "RandomForest": {"data": {}, "score": 0.95, "status": "MSS"}
                }
            }
        },
        "name": "splA",
        "results": {
            "RandomForest": {"score": 0.95, "status": "MSS"}
        }
    },
    {
        "loci": {
            "11:102322740-102322841": {
                "name": "NR27",
                "position": "11:102322740-102322841",
                "results": {
                    "RandomForest": {"data": {}, "score": 1.0, "status": "MSI"}
                }
            },
            "2:47414383-47414484": {
                "name": "BAT26",
                "position": "2:47414383-47414484",
                "results": {
                    "RandomForest": {"data": {}, "score": 0.95, "status": "MSI"}
                }
            }
        },
        "name": "splB",
        "results": {
            "RandomForest": {"score": 0.98, "status": "MSI"}
        }
    }
]""")

    def tearDown(self):
        for curr_file in [self.tmp_in_msings, self.tmp_in_msisensor, self.tmp_in_random, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test(self):
        cmd = [
            "microsatMergeResults.py",
            "--inputs-reports", self.tmp_in_msings, self.tmp_in_msisensor, self.tmp_in_random,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        expect = [
            {
                "loci": {
                    "11:102322740-102322841": {
                        "name": "NR27",
                        "position": "11:102322740-102322841",
                        "results": {
                            "mSINGS": {"data": {}, "score": None, "status": "MSS"},
                            "MSISensor": {"data": {}, "score": 1.0, "status": "MSS"},
                            "RandomForest": {"data": {}, "score": 0.99, "status": "MSS"}
                        }
                    },
                    "2:47414383-47414484": {
                        "name": "BAT26",
                        "position": "2:47414383-47414484",
                        "results": {
                            "mSINGS": {"data": {}, "score": None, "status": "MSS"},
                            "MSISensor": {"data": {}, "score": 1.0, "status": "MSS"},
                            "RandomForest": {"data": {}, "score": 0.95, "status": "MSS"}
                        }
                    }
                },
                "name": "splA",
                "results": {
                    "mSINGS": {"method": None, "param": None, "version": None, "score": 1.0, "status": "MSS"},
                    "MSISensor": {"method": None, "param": None, "version": None, "score": 1.0, "status": "MSS"},
                    "RandomForest": {"method": None, "param": None, "version": None, "score": 0.95, "status": "MSS"}
                }
            },
            {
                "loci": {
                    "11:102322740-102322841": {
                        "name": "NR27",
                        "position": "11:102322740-102322841",
                        "results": {
                            "mSINGS": {"data": {}, "score": None, "status": "MSS"},
                            "MSISensor": {"data": {}, "score": 0.9, "status": "MSI"},
                            "RandomForest": {"data": {}, "score": 1.0, "status": "MSI"}
                        }
                    },
                    "2:47414383-47414484": {
                        "name": "BAT26",
                        "position": "2:47414383-47414484",
                        "results": {
                            "mSINGS": {"data": {}, "score": None, "status": "MSI"},
                            "MSISensor": {"data": {}, "score": 0.8, "status": "MSI"},
                            "RandomForest": {"data": {}, "score": 0.95, "status": "MSI"}
                        }
                    }
                },
                "name": "splB",
                "results": {
                    "mSINGS": {"method": None, "param": None, "version": None, "score": 0.5, "status": "MSI"},
                    "MSISensor": {"method": None, "param": None, "version": None, "score": 0.95, "status": "MSI"},
                    "RandomForest": {"method": None, "param": None, "version": None, "score": 0.98, "status": "MSI"}
                }
            }
        ]
        with open(self.tmp_out) as reader:
            obs = json.load(reader)
        for elt_obs, elt_expect in zip(obs, expect):
            self.assertDictEqual(elt_obs, elt_expect)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
