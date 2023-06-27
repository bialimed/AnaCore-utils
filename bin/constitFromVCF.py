#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import argparse
from anacore.annotVcf import AnnotVCFIO, getAlleleRecord
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def extendCluster(AF_by_spl, clstr_samples, extend_rate=0.1):
    """
    Extend current cluster with superior neighbors AF. In each iteration, cluster
    is extended from is maximum allele frequency to the following if
    AF_neighbor <= clutr_max_AF + clutr_max_AF * extend_rate.

    :param AF_by_spl: Allele frequency by sample.
    :type AF_by_spl: dict
    :param clstr_samples: Samples in seed cluster.
    :type clstr_samples: list
    :param extend_rate: Maximum allele frequency increasing rate to aggregate the following frequency.
    :type extend_rate: float
    :return: Extension trace.
    :rtype: dict
    """
    clstr_max_AF = max([AF_by_spl[spl] for spl in clstr_samples])
    extended_threshold = clstr_max_AF + clstr_max_AF * extend_rate
    asc_observation = [{"spl": spl, "af": af} for spl, af in sorted(AF_by_spl.items(), key=lambda item: item[1]) if af > clstr_max_AF]
    trace = {"from": clstr_max_AF, "to": None, "ini_count": len(clstr_samples), "count": None}
    for obs in asc_observation:
        if obs["af"] >= clstr_max_AF and obs["af"] <= extended_threshold:
            clstr_samples.append(obs["spl"])
            clstr_max_AF = obs["af"]
            extended_threshold = clstr_max_AF + clstr_max_AF * extend_rate
            trace["to"] = obs["af"]
        else:
            break
    if trace["to"]:
        trace["count"] = len(clstr_samples)
    return trace


def getAFBySpl(reader, curr_allele, args):
    nb_usable_spl = 0  # Nb samples with valid DP at the variant position
    nb_support_spl = 0  # Nb samples containing the variant
    AF_by_spl = dict()
    for curr_spl in reader.samples:
        curr_DP = curr_allele.getDP(curr_spl)
        if curr_DP is None:
            raise Exception('The DP for the variant "{}" in sample "{}" must be set when you merge samples.'.format(curr_allele.getName(), curr_spl))
        elif curr_DP > args.min_DP:  # Skip samples with limited DP at the variant position
            nb_usable_spl += 1
            curr_AF = curr_allele.getAltAF(curr_spl)[0]
            if curr_AF is not None and curr_AF > 0:  # The sample contain the variant
                nb_support_spl += 1
                AF_by_spl[curr_spl] = curr_AF
    return AF_by_spl, nb_usable_spl, nb_support_spl


def getBiggerCluster(freq_by_spl, max_freq_diff):
    """
    Return the bigger group of samples with similar frequencies.

    :param freq_by_spl: The frequency by sample.
    :type freq_by_spl: dict
    :param max_freq_diff: The maximum difference of frequencies in a cluster.
    :type max_freq_diff: float
    :return: The list of samples in bigger cluster. The samples are sorted by frequency.
    :rtype: list
    """
    bigger_cluster_size = 0
    bigger_clster_samples = list()
    asc_freq_spl = sorted(freq_by_spl, key=freq_by_spl.get)
    spl_len = len(asc_freq_spl)
    for start_idx, start_spl in enumerate(asc_freq_spl):
        clstr_size = 0
        clstr_samples = list()  # The start sample itself will be added in comparison below
        start_AF = freq_by_spl[start_spl]
        if len(asc_freq_spl) - start_idx > bigger_cluster_size:
            curr_idx = start_idx
            out_of_range = False
            while curr_idx < spl_len - 1 and not out_of_range:
                curr_spl = asc_freq_spl[curr_idx]
                curr_spl_AF = freq_by_spl[curr_spl]
                if curr_spl_AF - start_AF <= max_freq_diff:
                    clstr_size += 1
                    clstr_samples.append(curr_spl)
                else:
                    out_of_range = True
                curr_idx += 1
        if clstr_size > bigger_cluster_size:
            bigger_cluster_size = clstr_size
            bigger_clster_samples = clstr_samples
    return bigger_clster_samples


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Find constitutive variants in population of samples coming from one VCF. An constitutive variant is a variant coming from polymerase error rate at this position, or workflow artifact, ... . It is common in almost all samples with an allele frequency very similar in each samples. If this variant is present in a superior frequency in a particular sample that can signify that this variant is not an artefact in this sample only.')
    parser.add_argument('-n', '--noise-offset', type=float, default=0.01, help='Value added to the maximum constitutive frequency of the variant. [Default: %(default)s]')
    parser.add_argument('-e', '--extension-rate', type=float, help='Extention of cluster of constitutive frequencies with neighbor values. Cluster are iteratively extend with following upper frequency if AF_next <= clutr_max_AF + clutr_max_AF * extend_rate.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_thresholds = parser.add_argument_group('Thresolds')  # Thresholds
    group_thresholds.add_argument('-r', '--min-presence-ratio', type=float, default=0.75, help='Minimum ratio between samples with variant at the same AF and the samples with valid depth (see "--min-DP") for declare variant as constitutive. [Default: %(default)s]')
    group_thresholds.add_argument('-c', '--min-presence-count', type=int, default=4, help='Minimum number of samples with variant at the same AF for declare variant as constitutive. [Default: %(default)s]')
    group_thresholds.add_argument('-d', '--min-DP', type=int, default=150, help='Minimum number of reads at the variant position for use the sample in constitutive variant evaluation. [Default: %(default)s]')
    group_thresholds.add_argument('-a', '--max-AF-var', type=float, default=0.015, help='Maximum variation between samples AF in constitutive variant. Constitutive variant have an AF similar between all samples, for determine if a variant is constitutive only the samples where the variant has an AF similar in comparison to this value are counted as support. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file containing the constitutive variants (format: TSV).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with AnnotVCFIO(args.input_variants) as reader:
        with open(args.output_variants, "w") as writer:
            # Header
            writer.write(
                "## PARAMETERS: {}\n".format(
                    " ".join(sys.argv)
                )
            )
            writer.write("## VERSION: {}\n".format(__version__))
            writer.write("\t".join([
                "#Chromosome",
                "Position",
                "Reference_allele",
                "Alternative_allele",
                "Noise_rate",
                "Nb_input_spl",
                "Nb_usable_spl",
                "Nb_support_spl",
                "Nb_constit_spl",
                "%_support_spl",
                "%_constit_spl",
                "%_under_threshold_spl",
                "filters",
                "Gene",
                "HGVSc",
                "HGVSp",
                "Constit_spl",
                "Constit_AF",
                "AF_over_noise_rate"
            ]) + "\n")
            # Records
            nb_spl = len(reader.samples)
            for record in reader:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(reader, record, idx)
                    AF_by_spl, nb_usable_spl, nb_support_spl = getAFBySpl(reader, curr_allele, args)
                    if nb_usable_spl > 0:
                        if nb_support_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_spl >= args.min_presence_count:  # Prefilter by support number for reduce processing time
                            clstr_samples = getBiggerCluster(AF_by_spl, args.max_AF_var)
                            if args.extension_rate:
                                extension_log = extendCluster(AF_by_spl, clstr_samples, args.extension_rate)
                                if extension_log["to"]:
                                    log.debug(
                                        "Extend {} cluster from {} to {} ({} samples to {})".format(
                                            curr_allele.getName(),
                                            extension_log["from"], extension_log["to"],
                                            extension_log["ini_count"], extension_log["count"]
                                        )
                                    )
                            nb_support_constit_spl = len(clstr_samples)
                            if nb_support_constit_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_constit_spl >= args.min_presence_count:  # The variant is constitutive
                                clstr_AF = [AF_by_spl[spl] for spl in sorted(clstr_samples)]
                                noise_threshold = max(clstr_AF)
                                if args.extension_rate:
                                    noise_threshold += noise_threshold * args.extension_rate
                                noise_threshold += args.noise_offset
                                noise_threshold = min(noise_threshold, 1)
                                AF_over_threshold = [AF for spl, AF in AF_by_spl.items() if AF > noise_threshold]
                                nb_under_threshold = nb_support_spl - len(AF_over_threshold)
                                genes = set()
                                hgvs_c = set()
                                hgvs_p = set()
                                if "ANN" in record.info:
                                    for annot in record.info["ANN"]:
                                        genes.add(annot["SYMBOL"])
                                        hgvs_c.add(annot["HGVSc"])
                                        hgvs_p.add(annot["HGVSp"])
                                genes = sorted(["" if elt is None else elt for elt in genes])
                                hgvs_c = sorted(["" if elt is None else elt for elt in hgvs_c])
                                hgvs_p = sorted(["" if elt is None else elt for elt in hgvs_p])
                                writer.write("{}\t{}\t{}\t{}\t{:.5f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    curr_allele.chrom,
                                    curr_allele.pos,
                                    curr_allele.ref,
                                    curr_allele.alt[0],
                                    noise_threshold,
                                    nb_spl,
                                    nb_usable_spl,
                                    nb_support_spl,
                                    nb_support_constit_spl,
                                    nb_support_spl / nb_usable_spl,
                                    nb_support_constit_spl / nb_support_spl,
                                    nb_under_threshold / nb_support_spl,
                                    ";".join(sorted(list(set(record.filter)))),
                                    ";".join(genes),
                                    ";".join(hgvs_c),
                                    ";".join(hgvs_p),
                                    ";".join(sorted(clstr_samples)),
                                    ";".join(map(str, sorted(clstr_AF))),
                                    ";".join(map(str, sorted(AF_over_threshold)))
                                ))
    log.info("End of job")
