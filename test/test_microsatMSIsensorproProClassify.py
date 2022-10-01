#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2022 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from copy import deepcopy
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
class TestMSIsensorproProClassify(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())
        self.tmp_in = os.path.join(tmp_folder, unique_id + "_in.json")
        with open(self.tmp_in, "w") as writer:
            writer.write("""[
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 6, "14": 3, "15": 6, "16": 16, "17": 44, "18": 56, "19": 94, "20": 127, "21": 121, "22": 73, "23": 67, "24": 35, "25": 15, "26": 62, "27": 1}, "mode": "reads"}}, "score": null, "status": null}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 8, "14": 13, "15": 26, "16": 72, "17": 181, "18": 349, "19": 463, "20": 472, "21": 602, "22": 62, "23": 3}, "mode": "reads"}}, "score": null, "status": null}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"11": 6, "12": 27, "13": 113, "14": 425, "15": 1019, "16": 1761, "17": 669, "18": 64, "19": 4, "20": 2, "21": 5}, "mode": "reads"}}, "score": null, "status": null}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"12": 2, "14": 1, "15": 12, "16": 21, "17": 35, "18": 78, "19": 190, "20": 331, "21": 705, "22": 161, "23": 79, "24": 45, "25": 11, "26": 5}, "mode": "reads"}}, "score": null, "status": null}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"11": 6, "12": 1, "13": 1, "14": 15, "15": 23, "16": 52, "17": 74, "18": 96, "19": 145, "20": 186, "21": 208, "22": 190, "23": 124, "24": 95, "25": 50, "26": 37, "27": 54, "28": 3, "29": 1}, "mode": "reads"}}, "score": null, "status": null}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"14": 3, "15": 8, "16": 6, "17": 10, "18": 46, "19": 95, "20": 170, "21": 145, "22": 82, "23": 154, "24": 8, "25": 2}, "mode": "reads"}}, "score": null, "status": null}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"14": 1, "15": 6, "16": 8, "17": 34, "18": 41, "19": 95, "20": 187, "21": 235, "22": 230, "23": 176, "24": 161, "25": 1203, "26": 62, "27": 15, "28": 1}, "mode": "reads"}}, "score": null, "status": null}}}
        },
        "name": "21T054986",
        "results": {}
    },
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"19": 4, "20": 6, "21": 3, "22": 6, "23": 11, "24": 4, "25": 2, "26": 22, "27": 1, "28": 1}, "mode": "reads"}}, "score": null, "status": null}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"17": 8, "18": 46, "19": 141, "20": 246, "21": 540, "22": 154, "23": 29, "24": 6}, "mode": "reads"}}, "score": null, "status": null}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"10": 3, "11": 12, "12": 52, "13": 181, "14": 642, "15": 1675, "16": 2756, "17": 788, "18": 65, "19": 3}, "mode": "reads"}}, "score": null, "status": null}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 2, "16": 1, "17": 4, "18": 11, "19": 10, "20": 25, "21": 107, "22": 36, "23": 34, "24": 10, "26": 1}, "mode": "reads"}}, "score": null, "status": null}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"17": 1, "18": 6, "19": 8, "20": 17, "21": 22, "22": 35, "23": 41, "24": 54, "25": 54, "26": 21, "27": 29}, "mode": "reads"}}, "score": null, "status": null}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"18": 2, "20": 1, "21": 4, "22": 4, "23": 17, "24": 2, "25": 2, "26": 1}, "mode": "reads"}}, "score": null, "status": null}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"18": 1, "19": 9, "20": 24, "21": 42, "22": 42, "23": 61, "24": 101, "25": 526, "26": 12, "27": 1, "28": 2}, "mode": "reads"}}, "score": null, "status": null}}}
        },
        "name": "21T017962",
        "results": {}
    }
]""")
        self.tmp_model = os.path.join(tmp_folder, unique_id + "_model.json")
        with open(self.tmp_model, "w") as writer:
            writer.write("""[
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.4275316740048522, "pro_q": 8.985533291400844e-05}, "lengths": {"ct_by_len": {"9": 4, "10": 7, "11": 44, "12": 99, "13": 43, "14": 67, "15": 73, "16": 13, "19": 2, "20": 7, "21": 9, "22": 15, "23": 17, "24": 8, "25": 8, "26": 11, "27": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 14, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.11407256264745265, "pro_q": 0.002423005802461264}, "lengths": {"ct_by_len": {"11": 1, "12": 2, "13": 18, "14": 27, "15": 56, "16": 49, "17": 88, "18": 57, "19": 95, "20": 134, "21": 185, "22": 28, "23": 5}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.13108990457056754, "pro_q": 0.003631727388633466}, "lengths": {"ct_by_len": {"9": 2, "10": 3, "11": 42, "12": 157, "13": 352, "14": 150, "15": 124, "16": 213, "17": 398, "18": 62, "19": 11, "20": 2, "21": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 8, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.10785241248817408, "pro_q": 0.024503311258278145}, "lengths": {"ct_by_len": {"12": 1, "13": 7, "14": 13, "15": 36, "16": 46, "17": 50, "18": 57, "19": 65, "20": 37, "21": 64, "22": 27, "23": 42, "24": 37, "25": 8, "26": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 13, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.4874201006391949, "pro_q": 0.0}, "lengths": {"ct_by_len": {"9": 8, "10": 28, "11": 94, "12": 196, "13": 133, "14": 123, "15": 138, "16": 31, "17": 1, "18": 2, "19": 4, "20": 1, "21": 6, "22": 9, "23": 4, "24": 6, "25": 19, "26": 10, "27": 4}, "mode": "reads"}, "mSINGS": {"nb_peaks": 9, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.14095500459136823, "pro_q": 0.002066115702479339}, "lengths": {"ct_by_len": {"12": 1, "13": 1, "14": 2, "15": 18, "16": 26, "17": 10, "18": 15, "19": 6, "20": 9, "21": 19, "22": 25, "23": 49, "24": 7, "25": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.14271930529870547, "pro_q": 0.004985125030151966}, "lengths": {"ct_by_len": {"13": 5, "14": 8, "15": 28, "16": 53, "17": 41, "18": 33, "19": 21, "20": 11, "21": 8, "22": 13, "23": 16, "24": 27, "25": 191, "26": 21, "27": 16, "28": 3}, "mode": "reads"}, "mSINGS": {"nb_peaks": 12, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}}
        },
        "name": "22T004446",
        "results": {}
    },
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.2428755419304515, "pro_q": 0.00012697952001741432}, "lengths": {"ct_by_len": {"12": 3, "13": 16, "14": 44, "15": 96, "16": 136, "17": 251, "18": 333, "19": 257, "20": 218, "21": 152, "22": 148, "23": 158, "24": 133, "25": 73, "26": 96, "27": 5, "28": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 13, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.10520380715911443, "pro_q": 0.005183116076970825}, "lengths": {"ct_by_len": {"11": 2, "12": 5, "13": 14, "14": 100, "15": 253, "16": 423, "17": 578, "18": 527, "19": 520, "20": 726, "21": 1018, "22": 341, "23": 58, "24": 12, "25": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.11974405850091407, "pro_q": 0.00047178156513534233}, "lengths": {"ct_by_len": {"10": 13, "11": 87, "12": 265, "13": 865, "14": 1440, "15": 1954, "16": 2618, "17": 675, "18": 56, "19": 1, "20": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 6, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.22308674701381212, "pro_q": 0.0151214329408987}, "lengths": {"ct_by_len": {"9": 2, "10": 39, "11": 175, "12": 280, "13": 560, "14": 675, "15": 535, "16": 376, "17": 271, "18": 124, "19": 89, "20": 135, "21": 384, "22": 258, "23": 293, "24": 139, "25": 27, "26": 8}, "mode": "reads"}, "mSINGS": {"nb_peaks": 15, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.2966028795848434, "pro_q": 0.00015511965366011874}, "lengths": {"ct_by_len": {"9": 5, "10": 13, "11": 44, "12": 115, "13": 190, "14": 274, "15": 259, "16": 133, "17": 105, "18": 117, "19": 127, "20": 148, "21": 151, "22": 146, "23": 179, "24": 175, "25": 172, "26": 124, "27": 141, "28": 5, "29": 3}, "mode": "reads"}, "mSINGS": {"nb_peaks": 17, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.1886397347526394, "pro_q": 0.0026175726376406945}, "lengths": {"ct_by_len": {"11": 2, "12": 9, "13": 26, "14": 66, "15": 148, "16": 126, "17": 91, "18": 49, "19": 47, "20": 40, "21": 57, "22": 68, "23": 219, "24": 35, "25": 8, "26": 3}, "mode": "reads"}, "mSINGS": {"nb_peaks": 12, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.11100356597045338, "pro_q": 0.0023773136355917813}, "lengths": {"ct_by_len": {"9": 4, "10": 2, "12": 2, "13": 9, "14": 16, "15": 64, "16": 158, "17": 292, "18": 397, "19": 420, "20": 265, "21": 135, "22": 129, "23": 205, "24": 311, "25": 2078, "26": 158, "27": 48, "28": 2, "29": 5}, "mode": "reads"}, "mSINGS": {"nb_peaks": 11, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSI"}}}
        },
        "name": "21T038725",
        "results": {}
    },
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.05369733205482543, "pro_q": 0.0012973094929212025}, "lengths": {"ct_by_len": {"15": 1, "16": 1, "17": 3, "18": 5, "19": 11, "20": 23, "21": 19, "22": 52, "23": 46, "24": 67, "25": 74, "26": 360, "27": 16, "28": 2, "29": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 7, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.03157356512846734, "pro_q": 0.009235587370573948}, "lengths": {"ct_by_len": {"14": 3, "15": 8, "16": 21, "17": 69, "18": 168, "19": 504, "20": 1002, "21": 1939, "22": 577, "23": 124, "24": 14}, "mode": "reads"}, "mSINGS": {"nb_peaks": 6, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.06374935087415613, "pro_q": 0.0026588194564652933}, "lengths": {"ct_by_len": {"10": 2, "11": 1, "12": 32, "13": 174, "14": 499, "15": 1665, "16": 3504, "17": 2270, "18": 277, "19": 44, "20": 3, "21": 1, "23": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 5, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.011713799471340117, "pro_q": 0.04797769489807003}, "lengths": {"ct_by_len": {"14": 2, "16": 6, "17": 13, "18": 52, "19": 99, "20": 197, "21": 676, "22": 590, "23": 615, "24": 199, "25": 46, "26": 7, "27": 1, "29": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 8, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.12249916392894379, "pro_q": 0.00037377294277339524}, "lengths": {"ct_by_len": {"12": 1, "13": 1, "14": 1, "16": 5, "17": 7, "18": 25, "19": 80, "20": 82, "21": 123, "22": 185, "23": 276, "24": 330, "25": 320, "26": 224, "27": 210, "28": 7, "29": 3, "30": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.045705010465165945, "pro_q": 0.00478407586177438}, "lengths": {"ct_by_len": {"13": 1, "14": 1, "15": 4, "16": 4, "17": 12, "18": 23, "19": 56, "20": 67, "21": 116, "22": 147, "23": 489, "24": 78, "25": 13, "27": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 6, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.014747524658393355, "pro_q": 0.005112221821015223}, "lengths": {"ct_by_len": {"13": 1, "14": 1, "16": 2, "17": 3, "18": 7, "19": 12, "20": 20, "21": 68, "22": 95, "23": 156, "24": 397, "25": 3044, "26": 260, "27": 97, "28": 18, "29": 5, "30": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 4, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}}
        },
        "name": "21T007948",
        "results": {}
    },
    {
        "loci": {
            "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.10686353829557713, "pro_q": 0.0008090614886731392}, "lengths": {"ct_by_len": {"15": 2, "16": 8, "17": 1, "18": 7, "19": 17, "20": 44, "21": 51, "22": 75, "23": 88, "24": 81, "25": 54, "26": 136, "27": 2, "28": 2, "29": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.025278012117493674, "pro_q": 0.011764705882352941}, "lengths": {"ct_by_len": {"14": 1, "15": 2, "16": 16, "17": 31, "18": 115, "19": 259, "20": 562, "21": 1431, "22": 543, "23": 100, "24": 8}, "mode": "reads"}, "mSINGS": {"nb_peaks": 6, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.03269050233105797, "pro_q": 0.004701895078524035}, "lengths": {"ct_by_len": {"10": 2, "12": 9, "13": 55, "14": 209, "15": 628, "16": 1947, "17": 3979, "18": 472, "19": 55, "20": 3}, "mode": "reads"}, "mSINGS": {"nb_peaks": 5, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.009708090816190519, "pro_q": 0.052494779401963836}, "lengths": {"ct_by_len": {"11": 1, "13": 3, "15": 3, "16": 9, "17": 11, "18": 25, "19": 53, "20": 115, "21": 621, "22": 402, "23": 490, "24": 225, "25": 61, "26": 10, "27": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 7, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.15256566202639996, "pro_q": 0.00013515339910798757}, "lengths": {"ct_by_len": {"12": 2, "14": 4, "15": 2, "16": 3, "17": 14, "18": 54, "19": 58, "20": 109, "21": 158, "22": 233, "23": 303, "24": 312, "25": 222, "26": 76, "27": 90, "28": 2, "29": 2}, "mode": "reads"}, "mSINGS": {"nb_peaks": 10, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.039659079712286924, "pro_q": 0.003990540939994088}, "lengths": {"ct_by_len": {"15": 1, "16": 2, "17": 10, "18": 12, "19": 29, "20": 64, "21": 110, "22": 135, "23": 448, "24": 58, "25": 7, "26": 3}, "mode": "reads"}, "mSINGS": {"nb_peaks": 6, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}},
            "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"model": {"data": {"MSIsensor-pro": {"pro_p": 0.015085387276118727, "pro_q": 0.002381426682844673}, "lengths": {"ct_by_len": {"13": 1, "17": 5, "18": 7, "19": 16, "20": 37, "21": 66, "22": 100, "23": 159, "24": 402, "25": 3400, "26": 175, "27": 30, "28": 8, "29": 1}, "mode": "reads"}, "mSINGS": {"nb_peaks": 3, "peak_height_cutoff": 0.05}}, "score": null, "status": "MSS"}}}
        },
        "name": "21T000980",
        "results": {}
    }
]""")
        self.tmp_out = os.path.join(tmp_folder, unique_id + "_out.json")
        self.expected = [
            {
                "loci": {
                    "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 6, "14": 3, "15": 6, "16": 16, "17": 44, "18": 56, "19": 94, "20": 127, "21": 121, "22": 73, "23": 67, "24": 35, "25": 15, "26": 62, "27": 1}, "mode": "reads"}, "pro_p": 0.20061450442337236, "pro_q": 5.297451925623775e-05}, "score": 0.999643, "status": "MSI"}}},
                    "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 8, "14": 13, "15": 26, "16": 72, "17": 181, "18": 349, "19": 463, "20": 472, "21": 602, "22": 62, "23": 3}, "mode": "reads"}, "pro_p": 0.08111704936732926, "pro_q": 0.0014364477492131223}, "score": 1.0, "status": "MSI"}}},
                    "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"11": 6, "12": 27, "13": 113, "14": 425, "15": 1019, "16": 1761, "17": 669, "18": 64, "19": 4, "20": 2, "21": 5}, "mode": "reads"}, "pro_p": 0.08172076944042, "pro_q": 0.0014057636308866352}, "score": 1.0, "status": "MSS"}}},
                    "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"12": 2, "14": 1, "15": 12, "16": 21, "17": 35, "18": 78, "19": 190, "20": 331, "21": 705, "22": 161, "23": 79, "24": 45, "25": 11, "26": 5}, "mode": "reads"}, "pro_p": 0.036031243875808396, "pro_q": 0.014642067247123379}, "score": 1.0, "status": "MSI"}}},
                    "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"11": 6, "12": 1, "13": 1, "14": 15, "15": 23, "16": 52, "17": 74, "18": 96, "19": 145, "20": 186, "21": 208, "22": 190, "23": 124, "24": 95, "25": 50, "26": 37, "27": 54, "28": 3, "29": 1}, "mode": "reads"}, "pro_p": 0.2271985198084458, "pro_q": 0.00013604701784936875}, "score": 0.999999, "status": "MSI"}}},
                    "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"14": 3, "15": 8, "16": 6, "17": 10, "18": 46, "19": 95, "20": 170, "21": 145, "22": 82, "23": 154, "24": 8, "25": 2}, "mode": "reads"}, "pro_p": 0.1004231479825973, "pro_q": 0.0007151796888968353}, "score": 1.0, "status": "MSI"}}},
                    "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"14": 1, "15": 6, "16": 8, "17": 34, "18": 41, "19": 95, "20": 187, "21": 235, "22": 230, "23": 176, "24": 161, "25": 1203, "26": 62, "27": 15, "28": 1}, "mode": "reads"}, "pro_p": 0.07076622742801367, "pro_q": 0.001545469334634781}, "score": 1.0, "status": "MSI"}}}
                },
                "name": "21T054986",
                "results": {
                    "MSIsensor-pro_pro": {"method": "MSIsensor-pro_pro", "param": {"aggregation_method": "instability ratio", "instability_threshold": 0.2, "min_voting_loci": 0.5}, "score": 0.85709, "status": "MSI", "version": "1.0.0"}
                }
            },
            {
                "loci": {
                    "11:102322777-102322803": {"name": "NR27", "position": "11:102322777-102322803", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"19": 4, "20": 6, "21": 3, "22": 6, "23": 11, "24": 4, "25": 2, "26": 22, "27": 1, "28": 1}, "mode": "reads"}, "pro_p": 0.09341010876519514, "pro_q": 0.0019193857965451055}, "score": 0.989567, "status": "MSS"}}},
                    "11:125620870-125620891": {"name": "NR22", "position": "11:125620870-125620891", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"17": 8, "18": 46, "19": 141, "20": 246, "21": 540, "22": 154, "23": 29, "24": 6}, "mode": "reads"}, "pro_p": 0.02814516129032258, "pro_q": 0.009274193548387096}, "score": 1.0, "status": "MSS"}}},
                    "13:31148483-31148500": {"name": "HT17", "position": "13:31148483-31148500", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"10": 3, "11": 12, "12": 52, "13": 181, "14": 642, "15": 1675, "16": 2756, "17": 788, "18": 65, "19": 3}, "mode": "reads"}, "pro_p": 0.0866863342215455, "pro_q": 0.0006756756756756757}, "score": 1.0, "status": "MSS"}}},
                    "14:23183137-23183158": {"name": "NR21", "position": "14:23183137-23183158", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"13": 2, "16": 1, "17": 4, "18": 11, "19": 10, "20": 25, "21": 107, "22": 36, "23": 34, "24": 10, "26": 1}, "mode": "reads"}, "pro_p": 0.022115384615384617, "pro_q": 0.02673076923076923}, "score": 1.0, "status": "MSI"}}},
                    "2:47414420-47414447": {"name": "BAT26", "position": "2:47414420-47414447", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"17": 1, "18": 6, "19": 8, "20": 17, "21": 22, "22": 35, "23": 41, "24": 54, "25": 54, "26": 21, "27": 29}, "mode": "reads"}, "pro_p": 0.1297582304526749, "pro_q": 0.0}, "score": 0.995897, "status": "MSS"}}},
                    "2:95183613-95183636": {"name": "NR24", "position": "2:95183613-95183636", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"18": 2, "20": 1, "21": 4, "22": 4, "23": 17, "24": 2, "25": 2, "26": 1}, "mode": "reads"}}, "score": None, "status": "Undetermined"}}},
                    "4:54732045-54732070": {"name": "BAT25", "position": "4:54732045-54732070", "results": {"MSIsensor-pro_pro": {"data": {"lengths": {"ct_by_len": {"18": 1, "19": 9, "20": 24, "21": 42, "22": 42, "23": 61, "24": 101, "25": 526, "26": 12, "27": 1, "28": 2}, "mode": "reads"}, "pro_p": 0.03397420296909224, "pro_q": 0.0009734728644439036}, "score": 1.0, "status": "MSI"}}}
                },
                "name": "21T017962",
                "results": {
                    "MSIsensor-pro_pro": {"method": "MSIsensor-pro_pro", "param": {"aggregation_method": "instability ratio", "instability_threshold": 0.2, "min_voting_loci": 0.5}, "score": 0.33576, "status": "MSI", "version": "1.0.0"}
                }
            }
        ]

    def tearDown(self):
        for curr_file in [self.tmp_in, self.tmp_model, self.tmp_out]:
            if os.path.exists(curr_file):
                os.remove(curr_file)

    def test_instabilityRatio(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--locus-weight-is-score",
            "--instability-ratio", "0.5",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        expected = deepcopy(self.expected)
        expected[0]["results"]["MSIsensor-pro_pro"]["param"]["instability_threshold"] = 0.5
        expected[1]["results"]["MSIsensor-pro_pro"] = {
            "method": "MSIsensor-pro_pro",
            "param": {"aggregation_method": "instability ratio", "instability_threshold": 0.5, "min_voting_loci": 0.5},
            "score": 0.66424,
            "status": "MSS",
            "version": "1.0.0"
        }
        for curr_obs, curr_expect in zip(observed, expected):
            self.assertDictEqual(curr_obs, curr_expect)

    def test_locusScore(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--locus-weight-is-score",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        for curr_obs, curr_expect in zip(observed, self.expected):
            self.assertDictEqual(curr_obs, curr_expect)

    def test_stdScore(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        expected = deepcopy(self.expected)
        expected[0]["results"]["MSIsensor-pro_pro"]["score"] = round(6 / 7, 5)
        expected[1]["results"]["MSIsensor-pro_pro"]["score"] = round(2 / 6, 5)
        for curr_obs, curr_expect in zip(observed, expected):
            self.assertDictEqual(curr_obs, curr_expect)

    def test_undeterminedScore(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--undetermined-weight", "0.3",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        expected = deepcopy(self.expected)
        expected[0]["results"]["MSIsensor-pro_pro"]["score"] = round(6 / 7, 5)
        expected[1]["results"]["MSIsensor-pro_pro"]["score"] = round(2 / (6 + 0.3), 5)
        for curr_obs, curr_expect in zip(observed, expected):
            self.assertDictEqual(curr_obs, curr_expect)

    def test_undeterminedAndLocusScore(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--locus-weight-is-score",
            "--undetermined-weight", "0.3",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        expected = deepcopy(self.expected)
        expected[1]["results"]["MSIsensor-pro_pro"]["score"] = round(
            (0.010433 + 0 + 0 + 1 + 0.004103 + 1) / 6.3,
            5
        )
        for curr_obs, curr_expect in zip(observed, expected):
            self.assertDictEqual(curr_obs, curr_expect)

    def test_votingLoci(self):
        cmd = [
            "microsatMSIsensorproProClassify.py",
            "--locus-weight-is-score",
            "--min-voting-loci", "0.9",
            "--input-evaluated", self.tmp_in,
            "--input-model", self.tmp_model,
            "--output-report", self.tmp_out
        ]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(self.tmp_out) as reader:
            observed = json.load(reader)
        expected = deepcopy(self.expected)
        expected[0]["results"]["MSIsensor-pro_pro"]["param"]["min_voting_loci"] = 0.9
        expected[1]["results"]["MSIsensor-pro_pro"] = {
            "method": "MSIsensor-pro_pro",
            "param": {"aggregation_method": "instability ratio", "instability_threshold": 0.2, "min_voting_loci": 0.9},
            "score": None,
            "status": "Undetermined",
            "version": "1.0.0"
        }
        for curr_obs, curr_expect in zip(observed, expected):
            self.assertDictEqual(curr_obs, curr_expect)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
