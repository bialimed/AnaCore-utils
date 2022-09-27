#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '3.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import json
import logging
import argparse
from anacore.msi.base import LocusClassifier, Status
from anacore.msi.locus import LocusRes
from anacore.msi.reportIO import ReportIO
from sklearn.tree import DecisionTreeClassifier as DecisionTree
from sklearn.neighbors import KNeighborsClassifier as KNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier as RandomForest
from sklearn.svm import SVC


########################################################################
#
# FUNCTIONS
#
########################################################################
class ClassifierParamsAction(argparse.Action):
    """Manages classifier-params parameters."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, json.loads(values))


class SklearnClassifier(LocusClassifier):
    def __init__(self, locus_id, method_name="MIAmS", model_method_name="model", clf="SVC", clf_params=None):
        if clf_params is None:
            clf_params = {}
        clf_obj = self._getClassifier(clf, clf_params)
        super().__init__(locus_id, method_name, clf_obj, model_method_name)

    def _getClassifier(self, clf, clf_params):
        clf_obj = None
        if clf == "SVC":  # The argument "probability" must be set to True to use predict_proba()
            clf_params["probability"] = True
            clf_params["gamma"] = "auto"
            clf_obj = SVC(**clf_params)
        elif clf == "KNeighbors":  # The KNeighbors does not accept the argument "random_state"
            if "n_neighbors" in clf_params:
                clf_params["n_neighbors"] = 2
            if "random_state" in clf_params:
                del clf_params["random_state"]
            clf_obj = KNeighbors(**clf_params)
        else:
            try:
                clf_obj = globals()[clf](**clf_params)
            except Exception:
                raise Exception('The classifier "{}" is not implemented in MIAmSClassifier.'.format(clf))
        return clf_obj


def process(args):
    """
    Predict classification (status and score) for all samples loci.

    :param args: The namespace extracted from the script arguments.
    :type args: Namespace
    """
    train_dataset = ReportIO.parse(args.input_model)
    test_dataset = ReportIO.parse(args.input_evaluated)
    # Classification by locus
    loci_ids = sorted(train_dataset[0].loci.keys())
    for locus_id in loci_ids:
        # Select the samples with a sufficient number of fragment to classify distribution
        evaluated_test_dataset = []
        for spl in test_dataset:
            locus = spl.loci[locus_id]
            locus_data = locus.results[args.data_method].data
            if args.data_method != args.status_method:  # Data come from another method
                locus_data = {"lengths": locus_data["lengths"]}
            locus.results[args.status_method] = LocusRes(
                Status.undetermined, None, locus_data
            )
            if locus_data["lengths"].getCount() >= args.min_depth:
                evaluated_test_dataset.append(spl)
        # Classify
        if len(evaluated_test_dataset) != 0:
            clf = SklearnClassifier(locus_id, args.status_method, "model", args.classifier, args.classifier_params)
            clf.fit(train_dataset)
            clf.set_status(evaluated_test_dataset)
    # Classification by sample
    for spl in test_dataset:
        spl.setStatusByInstabilityRatio(args.status_method, args.min_voting_loci, args.instability_ratio)
        spl.setScore(args.status_method, args.undetermined_weight, args.locus_weight_is_score)
    # Write output
    ReportIO.write(test_dataset, args.output_report)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Predict stability classes and scores for loci and samples using an sklearn classifer.')
    parser.add_argument('--data-method', help='The name of the method storing locus metrics and where the status will be set. [Default: classifier name]')
    parser.add_argument('--status-method', help='The name of the method storing locus metrics and where the status will be set. [Default: classifier name]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_locus = parser.add_argument_group('Locus classifier')  # Locus status
    group_locus.add_argument('-k', '--classifier', default="SVC", choices=["DecisionTree", "KNeighbors", "LogisticRegression", "RandomForest", "SVC"], help='The classifier used to predict loci status.')
    group_locus.add_argument('-p', '--classifier-params', action=ClassifierParamsAction, default={}, help='By default the classifier is used with these default parameters defined in scikit-learn. If you want change these parameters you use this option to provide them as json string. Example: {"n_estimators": 1000, "criterion": "entropy"} for RandmForest.')
    group_locus.add_argument('-f', '--min-depth', default=60, type=int, help='The minimum numbers of reads or fragments to determine the status. [Default: %(default)s]')
    group_locus.add_argument('-s', '--random-seed', default=None, type=int, help='The seed used by the random number generator in the classifier.')
    group_status = parser.add_argument_group('Sample consensus status')  # Sample status
    group_status.add_argument('-l', '--min-voting-loci', default=0.5, type=float, help='Minimum number of voting loci (stable + unstable) to determine the sample status. If the number of voting loci is lower than this value the status for the sample will be undetermined. [Default: %(default)s]')
    group_status.add_argument('-i', '--instability-ratio', default=0.2, type=float, help='If the ratio unstable/(stable + unstable) is superior than this value the status of the sample will be unstable otherwise it will be stable. [Default: %(default)s]')
    group_score = parser.add_argument_group('Sample prediction score')  # Sample score
    group_score.add_argument('-w', '--undetermined-weight', default=0, type=float, help='The weight of the undetermined loci in sample score calculation. [Default: %(default)s]')
    group_score.add_argument('-d', '--locus-weight-is-score', action='store_true', help='Use the prediction score of each locus as wheight of this locus in sample prediction score calculation. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-r', '--input-model', required=True, help='Path to the file containing the references samples used in learn step (format: MSIReport).')
    group_input.add_argument('-e', '--input-evaluated', required=True, help='Path to the file containing the samples with loci to classify (format: MSIReport).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-report', required=True, help='The path to the output file (format: MSIReport).')
    args = parser.parse_args()

    args.classifier_params["random_state"] = args.random_seed
    if args.data_method is None:
        args.data_method = args.classifier
    if args.status_method is None:
        args.status_method = args.classifier

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    process(args)
    log.info("End of job")
