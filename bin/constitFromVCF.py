#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2018 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '2.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import argparse
from anacore.annotVcf import AnnotVCFIO, getAlleleRecord, VCFIO
from anacore.sv import HashedSVIO
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
class Cluster:
    def __init__(self, af_by_spl, min_idx=None, max_idx=None):
        self.max_idx = max_idx
        self.min_idx = min_idx
        self.spl = sorted(af_by_spl, key=af_by_spl.get)
        self.af = [af_by_spl[spl] for spl in self.spl]

    def extend(self, extend_rate=0.1):
        """
        Extend current cluster with superior and inferior neighbors AF. In each
        iteration, cluster is extended from is maximum allele frequency to the
        following if AF_neighbor <= clutr_max_AF + clutr_max_AF * extend_rate.

        :param extend_rate: Maximum allele frequency increasing rate to aggregate the following frequency.
        :type extend_rate: float
        """
        # Extend max
        ini_max_af = self.max_af
        extended_threshold = ini_max_af + ini_max_af * extend_rate
        for curr_af in self.over_af:
            if curr_af <= extended_threshold:
                self.max_idx += 1
                extended_threshold = curr_af + curr_af * extend_rate
            else:
                break
        # Extend min
        ini_min_af = self.min_af
        extended_threshold = ini_min_af - ini_min_af * extend_rate
        for prev_af in self.under_af[::-1]:
            if prev_af >= extended_threshold:
                self.min_idx -= 1
                extended_threshold = prev_af - prev_af * extend_rate
            else:
                break

    @property
    def in_af(self):
        return self.af[self.min_idx:self.max_idx + 1]

    @property
    def in_spl(self):
        return self.spl[self.min_idx:self.max_idx + 1]

    @property
    def len(self):
        return self.max_idx - self.min_idx + 1

    @property
    def max_af(self):
        return self.af[self.max_idx]

    @property
    def min_af(self):
        return self.af[self.min_idx]

    @property
    def over_af(self):
        return self.af[self.max_idx + 1:]

    def setBiggerCluster(self, max_af_diff):
        """
        :param max_af_diff: Maximum difference of frequencies in a cluster.
        :type max_af_diff: float
        """
        bigger_clster = {
            "size": 0,
            "start": 0
        }
        pop_len = len(self.af)
        for start_idx, start_AF in enumerate(self.af):
            clstr_size = 1
            if pop_len - start_idx > bigger_clster["size"]:
                for end_AF in self.af[start_idx + 1:]:
                    if end_AF - start_AF <= max_af_diff:
                        clstr_size += 1
                    else:
                        break
            if clstr_size > bigger_clster["size"]:
                bigger_clster = {
                    "size": clstr_size,
                    "start": start_idx
                }
        self.min_idx = bigger_clster["start"]
        self.max_idx = bigger_clster["start"] + bigger_clster["size"] - 1

    @property
    def under_af(self):
        return self.af[:self.min_idx]


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


def getNoiseThreshold(clstr, extension_rate, noise_offset):
    noise_threshold = clstr.max_af
    if extension_rate:
        noise_threshold += noise_threshold * extension_rate
    noise_threshold += noise_offset
    noise_threshold = min(noise_threshold, 1)
    return noise_threshold


def getVarAnnot(record, annot_field):
    genes = set()
    hgvs_c = set()
    hgvs_p = set()
    if annot_field in record.info:
        for annot in record.info[annot_field]:
            genes.add(annot["SYMBOL"])
            hgvs_c.add(annot["HGVSc"])
            hgvs_p.add(annot["HGVSp"])
    return {
        "genes": sorted(["" if elt is None else elt for elt in genes]),
        "hgvs_c": sorted(["" if elt is None else elt for elt in hgvs_c]),
        "hgvs_p": sorted(["" if elt is None else elt for elt in hgvs_p])
    }


class LoggerAction(argparse.Action):
    """Manages logger level parameters (The value "INFO" becomes logging.info and so on)."""

    def __call__(self, parser, namespace, values, option_string=None):
        log_level = None
        if values == "DEBUG":
            log_level = logging.DEBUG
        elif values == "INFO":
            log_level = logging.INFO
        elif values == "WARNING":
            log_level = logging.WARNING
        elif values == "ERROR":
            log_level = logging.ERROR
        elif values == "CRITICAL":
            log_level = logging.CRITICAL
        setattr(namespace, self.dest, log_level)


def readerFactory(path, annot_field):
    with VCFIO(args.input_variants) as reader:
        if annot_field and annot_field in reader.info:
            final_reader = AnnotVCFIO(path, annot_field=annot_field)
        else:
            final_reader = VCFIO(args.input_variants)
    return final_reader


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Find constitutive variants in population of samples coming from one VCF. An constitutive variant is a variant coming from polymerase error rate at this position, or workflow artifact, ... . It is common in almost all samples with an allele frequency very similar in each samples. If this variant is present in a superior frequency in a particular sample that can signify that this variant is not an artefact in this sample only.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_thresholds = parser.add_argument_group('Thresolds')  # Thresholds
    group_thresholds.add_argument('-r', '--min-presence-ratio', type=float, default=0.75, help='Minimum ratio between samples with variant at the same AF and the samples with valid depth (see "--min-DP") for declare variant as constitutive. [Default: %(default)s]')
    group_thresholds.add_argument('-c', '--min-presence-count', type=int, default=4, help='Minimum number of samples with variant at the same AF for declare variant as constitutive. [Default: %(default)s]')
    group_thresholds.add_argument('-d', '--min-DP', type=int, default=150, help='Minimum number of reads at the variant position for use the sample in constitutive variant evaluation. [Default: %(default)s]')
    group_clstr = parser.add_argument_group('Cluster')  # Cluster
    group_clstr.add_argument('-e', '--extension-rate', type=float, help='Extention of cluster of constitutive frequencies with neighbor values. Cluster are iteratively extend with following upper frequency if AF_next <= clutr_max_AF + clutr_max_AF * extend_rate.')
    group_clstr.add_argument('-f', '--max-AF-var', type=float, default=0.015, help='Maximum variation between samples AF in constitutive variant. Constitutive variant have an AF similar between all samples, for determine if a variant is constitutive only the samples where the variant has an AF similar in comparison to this value are counted as support. [Default: %(default)s]')
    group_clstr.add_argument('-n', '--noise-offset', type=float, default=0.01, help='Value added to the maximum constitutive frequency of the variant. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file containing the constitutive variants (format: TSV).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with readerFactory(args.input_variants, args.annotation_field) as reader:
        with HashedSVIO(args.output_variants, "w", title_starter="#", metadata_starter="##") as writer:
            # Header
            writer.metadata.append("PARAMETERS: {}".format(" ".join(sys.argv)))
            writer.metadata.append("VERSION: {}".format(__version__))
            writer.titles = [
                "Chromosome", "Position", "Reference_allele", "Alternative_allele",
                "Noise_rate", "Nb_input_spl", "Nb_usable_spl", "Nb_support_spl",
                "Nb_constit_spl", "%_support_spl", "%_constit_spl", "%_under_threshold_spl",
                "Filters", "Gene", "HGVSc", "HGVSp", "Constit_spl", "Constit_AF",
                "AF_over_noise_rate", "All_AF"
            ]
            writer.writeHeader()
            # Records
            nb_spl = len(reader.samples)
            for record in reader:
                for idx in range(len(record.alt)):
                    curr_allele = getAlleleRecord(reader, record, idx)
                    AF_by_spl, nb_usable_spl, nb_support_spl = getAFBySpl(reader, curr_allele, args)
                    if nb_usable_spl > 0:
                        if nb_support_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_spl >= args.min_presence_count:  # Prefilter by support number for reduce processing time
                            clstr = Cluster(AF_by_spl)
                            clstr.setBiggerCluster(args.max_AF_var)
                            if clstr.len > 1:
                                log.debug("Biggest seed for {} contains {}/{} elements from {} to {}.".format(curr_allele.getName(), clstr.len, nb_usable_spl, clstr.min_af, clstr.max_af))
                            if args.extension_rate:
                                prev_nb_clstr_spl = clstr.len
                                clstr.extend(args.extension_rate)
                                if prev_nb_clstr_spl < clstr.len:
                                    log.debug("Biggest extended cluster for {} contains {}/{} elements from {} to {}.".format(curr_allele.getName(), clstr.len, nb_usable_spl, clstr.min_af, clstr.max_af))
                            nb_support_constit_spl = clstr.len
                            if nb_support_constit_spl / nb_usable_spl >= args.min_presence_ratio and nb_support_constit_spl >= args.min_presence_count:  # The variant is constitutive
                                noise_threshold = getNoiseThreshold(clstr, args.extension_rate, args.noise_offset)
                                AF_over_threshold = [af for af in clstr.over_af if af > noise_threshold]
                                nb_under_threshold = nb_support_spl - len(AF_over_threshold)
                                clstr_annot = getVarAnnot(record, args.annotation_field)
                                writer.write({
                                    "Chromosome": curr_allele.chrom,
                                    "Position": curr_allele.pos,
                                    "Reference_allele": curr_allele.ref,
                                    "Alternative_allele": curr_allele.alt[0],
                                    "Noise_rate": "{:.5f}".format(noise_threshold),
                                    "Nb_input_spl": nb_spl,
                                    "Nb_usable_spl": nb_usable_spl,
                                    "Nb_support_spl": nb_support_spl,
                                    "Nb_constit_spl": nb_support_constit_spl,
                                    "%_support_spl": "{:.5f}".format(nb_support_spl / nb_usable_spl),
                                    "%_constit_spl": "{:.5f}".format(nb_support_constit_spl / nb_support_spl),
                                    "%_under_threshold_spl": "{:.5f}".format(nb_under_threshold / nb_support_spl),
                                    "Filters": ";".join(sorted(list(set(record.filter)))),
                                    "Gene": ";".join(clstr_annot["genes"]),
                                    "HGVSc": ";".join(clstr_annot["hgvs_c"]),
                                    "HGVSp": ";".join(clstr_annot["hgvs_p"]),
                                    "Constit_spl": ";".join(sorted(clstr.in_spl)),
                                    "Constit_AF": ";".join(map(str, clstr.in_af)),
                                    "AF_over_noise_rate": ";".join(map(str, AF_over_threshold)),
                                    "All_AF": ";".join(map(str, clstr.af))
                                })
    log.info("End of job")
