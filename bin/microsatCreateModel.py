#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.bed import BEDIO
from anacore.msi.annot import addLociResToSpl, getLocusAnnotDict, MSIAnnot
from anacore.msi.base import Status
from anacore.msi.locus import Locus
from anacore.msi.msings import MSINGSEval
from anacore.msi.msisensorpro import ProEval
from anacore.msi.reportIO import ReportIO
import argparse
from copy import deepcopy
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def addMSINGSInfo(msi_samples, result_name, peak_height_cutoff):
    """
    Add {"nb_peaks": ..., "peak_height_cutoff": ...} from mSINGS in loci results.

    :param msi_samples: MSI samples.
    :type msi_samples: list of anacore.msi.sample.MSISample
    :param result_name: Name of the method used to store model data.
    :type result_name: str
    :param peak_height_cutoff: Minimum rate of the highest peak to consider a peak in size distribution.
    :type peak_height_cutoff: float
    """
    for curr_spl in msi_samples:
        for locus_id, locus in curr_spl.loci.items():
            if result_name in locus.results:
                locus_data = locus.results[result_name].data
                locus_data["mSINGS"] = {
                    "nb_peaks": MSINGSEval.getNbPeaks(locus_data["lengths"], peak_height_cutoff),
                    "peak_height_cutoff": peak_height_cutoff
                }


def addMSIsensorInfo(msi_samples, result_name):
    """
    Add {"pro_p": ..., "pro_q": ...}} from MSIsensor-pro in loci results.

    :param msi_samples: MSI samples.
    :type msi_samples: list of anacore.msi.sample.MSISample
    :param result_name: Name of the method used to store model data.
    :type result_name: str
    """
    for curr_spl in msi_samples:
        for locus_id, locus in curr_spl.loci.items():
            if result_name in locus.results:
                locus_data = locus.results[result_name].data
                pro_p, pro_q = ProEval.getSlippageScores(locus_data["lengths"], locus.end - locus.start + 1)
                locus_data["MSIsensor-pro"] = {
                    "pro_p": pro_p,
                    "pro_q": pro_q
                }


def getAggregatedSpl(in_reports):
    """
    Return one list of MSISample from several MSReport.

    :param in_reports: Pathes to the MSIReport files.
    :type in_reports: list
    :return: List of MSISample.
    :rtype: list
    """
    aggregated_spl = []
    for curr_report in in_reports:
        msi_samples = ReportIO.parse(curr_report)
        for curr_spl in msi_samples:
            aggregated_spl.append(curr_spl)
    return aggregated_spl


def populateLoci(msi_samples, ref_loci):
    """
    Add loci if they are missing in sample.

    :param msi_samples: The samples to populate.
    :type msi_samples: list of MSI samples
    :param ref_loci: The loci to add if they are missing in samples.
    :type ref_loci: str
    """
    for spl in msi_samples:
        for ref_locus in ref_loci:
            if ref_locus.position not in spl.loci:
                spl.addLocus(deepcopy(ref_locus))


def process(args):
    """
    Create training data for MSI classifiers. These references are stored in
    MSIReport format.

    :param args: The namespace extracted from the script arguments.
    :type args: Namespace
    """
    # Get method name from annotations file
    method_names = set()
    for record in MSIAnnot(args.input_loci_status):
        method_names.add(record["method_id"])
    if len(method_names) != 1:
        raise ValueError('The annotation file must contain only one value for method_id. The file "{}" contains {}.'.format(args.input_length_distributions, method_names))
    result_id = list(method_names)[0]
    # Get reference loci from targets file
    ref_loci = []
    with BEDIO(args.input_microsatellites) as FH_in:
        for record in FH_in:
            ref_loci.append(
                Locus(
                    "{}:{}-{}".format(record.chrom, record.start - 1, record.end),
                    record.name
                )
            )
    # Aggregate samples
    msi_samples = getAggregatedSpl(args.inputs_length_distributions)
    # Add locus result info
    data_by_spl = getLocusAnnotDict(args.input_loci_status)
    for curr_spl in msi_samples:
        addLociResToSpl(curr_spl, data_by_spl[curr_spl.name])
    # Filter locus results
    populateLoci(msi_samples, ref_loci)
    pruneResults(msi_samples, result_id, args.min_support)
    # Add classifiers data
    addMSINGSInfo(msi_samples, result_id, args.peak_height_cutoff)
    addMSIsensorInfo(msi_samples, result_id)
    # Display metrics
    writeStatusMetrics(msi_samples, result_id, args.output_info)
    # Write output
    ReportIO.write(msi_samples, args.output_model)


def pruneResults(msi_samples, result_id, min_support):
    """
    Remove LocusRes where the status is not determined (none or undetermined)
    and/or where the number of fragment used to determine status is lower than
    min_support.

    :param msi_samples: The pruned samples.
    :type msi_samples: list of MSISample
    :param result_id: The method on which the filters are applied.
    :type result_id: str
    :param min_support: The minimum number of support to keep a LocusRes in data.
    :type min_support: int
    """
    removed_spl_idx = list()
    for spl_idx, spl in enumerate(msi_samples):
        nb_results = 0
        for locus_id, msi_locus in spl.loci.items():
            if result_id in msi_locus.results:
                if msi_locus.results[result_id].status not in [Status.stable, Status.unstable]:
                    msi_locus.delResult(result_id)
                elif msi_locus.results[result_id].data["lengths"].getCount() < min_support:
                    msi_locus.delResult(result_id)
                else:
                    nb_results += 1
        if nb_results == 0:
            removed_spl_idx.append(spl_idx)
    for spl_idx in sorted(removed_spl_idx)[::-1]:
        del(msi_samples[spl_idx])


def writeStatusMetrics(msi_samples, result_id, out_summary):
    """
    Write the statistics of status by loci in population of samples.

    :param msi_samples: The samples processed.
    :type msi_samples: list of MSI samples
    :param result_id: Only the results of this methd are processed.
    :type result_id: str
    :param out_summary: Path to the output file.
    :type out_summary: str
    """
    status_by_locus = dict()
    locus_name_by_id = dict()
    authorized_status = sorted(Status.authorizedValues() - {Status.undetermined, Status.none})  # Removed pruned status
    # Get number of samples by status for each locus
    for spl in msi_samples:
        for locus_id, locus in spl.loci.items():
            locus_name_by_id[locus_id] = locus.name
            if locus_id not in status_by_locus:
                status_by_locus[locus_id] = {status: 0 for status in authorized_status}  # init each status for the current locus
            if result_id in locus.results:
                status = locus.results[result_id].status
                status_by_locus[locus_id][status] += 1
    # Write results
    nb_spl = len(msi_samples)
    report_status = authorized_status + ["Unused"]
    with open(out_summary, "w") as FH_out:
        FH_out.write("Nb retained samples: {}\n".format(nb_spl))
        print(
            "Locus_position", "Locus_name", "\t".join([str(status) for status in report_status]),
            sep="\t",
            file=FH_out
        )
        for locus_id, locus_name in locus_name_by_id.items():
            status_by_locus[locus_id]["Unused"] = nb_spl - sum(status_by_locus[locus_id].values())
            print(
                locus_id,
                locus_name,
                "\t".join(
                    [str(status_by_locus[locus_id][status]) for status in report_status]
                ),
                sep="\t",
                file=FH_out
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Create training data for MSI classifiers. These references are stored in MSIReport format. All the loci are represented in all samples but all the loci does not have a result (if the data does not fit filters criteria).')
    parser.add_argument('-p', '--peak-height-cutoff', default=0.05, type=float, help='[mSINGS] Minimum height to consider a peak in size distribution as rate of the highest peak. [Default: %(default)s]')
    parser.add_argument('-s', '--min-support', type=int, default=200, help='Minimum number of reads/fragments in size distribution to keep the result. The distribution must contains a sufficient amount of data to be representative of length distribution profile for the current locus. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-d', '--inputs-length-distributions', required=True, nargs='+', help='Path(es) to the file(s) evaluated in references creation process (format: MSIReport).')
    group_input.add_argument('-l', '--input-loci-status', required=True, help='Path to the file containing for each sample for each targeted locus the stability status (format: MSIAnnot). First line must be: sample<tab>locus_position<tab>method_id<tab>key<tab>value<tab>type. The method_id should be "model" and an example of line content is: H2291-1_S15<tab>4:55598140-55598290<tab>model<tab>status<tab>MSS<tab>str.')
    group_input.add_argument('-t', '--input-microsatellites', required=True, help='Path to file containing locations of the microsatellite of interest (format: BED).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-i', '--output-info', required=True, help='The path to the file describing the number of references by status for each locus (format: TSV).')
    group_output.add_argument('-m', '--output-model', required=True, help='The path to the file containing the references distribution for each locus (format: MSIReport).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    process(args)
    log.info("End of job")
