#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.2.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import re
import os
import sys
import logging
import argparse
import requests
import urllib.parse
from anacore.annotVcf import AnnotVCFIO
from anacore.sv import HashedSVIO
from anacore.hgvs import HGVS, RunMutalyzerDescription, RunMutalyzerLegend


########################################################################
#
# FUNCTIONS
#
########################################################################
def getHGVSgFromRec(record, acc_by_chrom=None, annotations_field="ANN"):
    """
    Return HGVSg from annotated VCF record.

    :param record: Annotated VCF record.
    :type record: anacore.annotVcf.VCFRecord
    :param acc_by_chrom: Chromosome RefSeq accession by chromosome name.
    :type acc_by_chrom: dict
    :param annotations_field: Field used to store annotations.
    :type annotations_field: str
    :return: HGVSg.
    :rtype: str
    """
    hgvs_g = set()
    for annot in record.info[annotations_field]:
        if "HGVSg" in annot and annot["HGVSg"] is not None and annot["HGVSg"] != "":
            hgvs_g.add(annot["HGVSg"])
    if len(hgvs_g) > 1:
        raise Exception("The variant {} is describes with several HGVSg: {}.".format(record.getName(), sorted(hgvs_g)))
    elif len(hgvs_g) != 0:
        hgvs_g = list(hgvs_g)[0]
        if not hgvs_g.startswith("NC"):
            chr_acc = acc_by_chrom[record.chrom]
            hgvs_g = hgvs_g.replace(record.chrom + ":", chr_acc + ":")
    else:  # Empty HGVSg
        hgvs_g = None
    return hgvs_g


def getHGVSByTr(res_data):
    """
    Return HGVSg, HGVSc/n and HGVSp by transcript base RefSeq accession from runMutalyzer[Light].

    :param res_data: Results from runMutalyzer[Light] (required fields: legend, genomicDescription, transcriptDescriptions and proteinDescriptions).
    :type res_data: list
    :return: HGVSg and HGVSg, HGVSc/n and HGVSp by transcript base RefSeq accession.
    :rtype: str, dict
    """
    legend = RunMutalyzerLegend(res_data["legend"])
    id_by_name = legend.getIdByName()
    prot_by_tr = legend.getProtBytr()
    new_HGVSg = res_data["genomicDescription"]
    HGVSc_by_tr = RunMutalyzerDescription(res_data["transcriptDescriptions"]).getByAccession(id_by_name)
    HGVSp_by_prot = RunMutalyzerDescription(res_data["proteinDescriptions"]).getByAccession(id_by_name)
    HGVS_by_tr = {}
    for tr_acc, HGVSc in HGVSc_by_tr.items():
        HGVSp = ""
        if tr_acc in prot_by_tr:  # if the mRNA is linked to a protein (for predicted mRNA and some variants out of CDS runMutalyzer does not return any protein)
            prot_acc = prot_by_tr[tr_acc]
            HGVSp = HGVSp_by_prot[prot_acc]
        tr_base_ac, tr_acc_version = tr_acc.split(".")
        tr_acc_version = int(tr_acc_version)
        if tr_base_ac not in HGVS_by_tr or tr_acc_version > HGVS_by_tr[tr_base_ac]["acc_version"]:
            HGVS_by_tr[tr_base_ac] = {
                "acc_version": tr_acc_version,
                "HGVSg": new_HGVSg,
                "HGVSc": HGVSc,
                "HGVSp": HGVSp
            }
    return new_HGVSg, HGVS_by_tr


def getConsistentHGVS(new_hgvs_str, old_hgvs_str):
    """
    Return the HGVS string coming from mutalyzer or the HGVS fixed on several known bug.

    :param new_hgvs_str: The HGVS coming from mutalyzer.
    :type new_hgvs_str: str
    :param old_hgvs_str: The previously HGVS.
    :type old_hgvs_str: str
    :return: The consensus HGVS string.
    :rtype: str
    """
    new_hgvs = None if new_hgvs_str == "" else HGVS.fromStr(new_hgvs_str)
    old_hgvs = None if old_hgvs_str == "" else HGVS.fromStr(old_hgvs_str)
    final_hgvs_str = new_hgvs_str
    if new_hgvs is None:
        final_hgvs_str = old_hgvs_str
    else:
        if new_hgvs.change == "?":
            final_hgvs_str = ""
        if old_hgvs_str != "":
            if new_hgvs.change == "?":
                final_hgvs_str = old_hgvs_str
            else:
                # No change because ref sequence from HGVS database contains the alt
                if "=" in new_hgvs.change and new_hgvs.type != "p":
                    match = re.search(r"\d([ATGCN]+>[ATGCN]+)$", old_hgvs)
                    if match is not None:
                        final_hgvs_str = final_hgvs_str.replace("=", match.group(1))
                # Keep consistency on accession version
                if str(old_hgvs) != str(new_hgvs):
                    if old_hgvs.change == new_hgvs.change and old_hgvs.accession.id == new_hgvs.accession.id:
                        final_hgvs_str = old_hgvs_str
    return final_hgvs_str


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


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Fix or add HGVSg, HGVSc and HGVSp on variants annotations. The HGVS used are based on mutalyzer.')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-m', '--mutalyzer-url', default="https://mutalyzer.nl", help='URL to the mutalizer server. [Default: %(default)s]')
    parser.add_argument('-p', '--proxy-url', help='URL to the proxy server if the http(s) connexions are only allowed through a proxy.')
    parser.add_argument('-s', '--assembly-version', default="GRCh38", help='Human genome assembly version used in alignment, variants calling and variants annotation. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF).')
    group_input.add_argument('-c', '--input-assembly-accessions', help='Path to the file describing link between chromosome name and RefSeq accession (format: TSV). The header must contain: sequence_id<tab>...<tab>RefSeq_accession where sequence_id is the name of the chromosome in CHROM column of the VCF. [Default: the chromosome name in VCF is the RefSeq accession]')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Get accession by chromosome ID
    acc_by_chrom = {}
    if args.input_assembly_accessions:
        with HashedSVIO(args.input_assembly_accessions, title_starter=None) as FH:
            for record in FH:
                acc_by_chrom[record["sequence_id"]] = record["RefSeq_accession"]

    # Write
    nb_records = {"analysed": 0, "fixed_HGVSg": 0, "fixed_HGVSc": 0, "fixed_HGVSp": 0, "contains_colloc_annot": 0}
    with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotations_field) as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotations_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            if "HGVSg" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSg")
            if "HGVSc" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSc")
            if "HGVSp" not in FH_out.ANN_titles:
                FH_out.ANN_titles.append("HGVSp")
            FH_out.writeHeader()
            # Records
            param_assembly = urllib.parse.quote(args.assembly_version, safe='')
            for record in FH_in:
                nb_records["analysed"] += 1
                if len(record.alt) > 1:
                    raise Exception("The record {} is multi-allelic.".format(record.getName()))
                # Get info from mutalyzer
                old_HGVSg = getHGVSgFromRec(record, acc_by_chrom, FH_in.annot_field)
                if old_HGVSg is None:
                    log.warning("The variant {} does not contain any HGVSg.".format(record.getName()))
                else:
                    param_hgvsg = urllib.parse.quote(old_HGVSg, safe='')
                    param_fields = urllib.parse.quote(",".join(["legend", "proteinDescriptions", "transcriptDescriptions", "genomicDescription"]), safe='')
                    url_request = '{}/json/runMutalyzerLight?build={};variant={};extra={}'.format(args.mutalyzer_url, param_assembly, param_hgvsg, param_fields)
                    log.debug(url_request)
                    response = requests.get(
                        url_request,
                        proxies=(None if args.proxy_url is None else {"https": args.proxy_url, "http": args.proxy_url})
                    )
                    if response.status_code != 200:
                        raise Exception("Request {} has failed.".format(url_request))
                    res_data = response.json()
                    if res_data["errors"] != 0:
                        log.warning("The variant {} cannot be standardized by mutalyzer because: {}".format(record.getName(), res_data["messages"]))
                    else:
                        # Store all HGVS by transcript base accession
                        mutalyzer_tr = {elt["id"].split(".")[0] for elt in res_data["legend"] if "id" in elt and not elt["id"].startswith("ENS")}
                        annot_tr = {annot["Feature"].split(".")[0] for annot in record.info[args.annotations_field] if annot["Feature"] and not annot["Feature"].startswith("ENS")}
                        if len(annot_tr - mutalyzer_tr) != 0:
                            log.warning("All the transcripts annotated for variant {} cannot be found in used version of mutalyzer. Missing transcripts: {}".format(record.getName(), sorted(annot_tr - mutalyzer_tr)))
                        new_HGVSg, HGVS_by_tr = getHGVSByTr(res_data)
                        # Update annotations
                        is_fixed_HGVSg = False
                        is_fixed_HGVSc = False
                        is_fixed_HGVSp = False
                        contains_colloc_annot = False
                        new_annot = []
                        for annot in record.info[args.annotations_field]:
                            if annot["Allele"] != record.alt[0]:  # Annotation come from a collocated alternative allele
                                contains_colloc_annot = True
                            else:  # Annotation come from the alternative allele
                                # HGVSg
                                old_HGVSg = "" if "HGVSg" not in annot or annot["HGVSg"] is None else annot["HGVSg"]
                                annot["HGVSg"] = getConsistentHGVS(new_HGVSg, old_HGVSg)
                                if old_HGVSg != annot["HGVSg"]:
                                    is_fixed_HGVSg = True
                                # HGVSc and p
                                if annot["Feature"]:
                                    tr_base_acc = annot["Feature"].split(".")[0]
                                    if tr_base_acc in HGVS_by_tr:
                                        # HGVSc
                                        old_HGVSc = "" if "HGVSc" not in annot or annot["HGVSc"] is None else annot["HGVSc"]
                                        annot["HGVSc"] = getConsistentHGVS(HGVS_by_tr[tr_base_acc]["HGVSc"], old_HGVSc)
                                        if old_HGVSc != annot["HGVSc"]:
                                            is_fixed_HGVSc = True
                                        # HGVSp
                                        old_HGVSp = "" if "HGVSp" not in annot or annot["HGVSp"] is None else annot["HGVSp"]
                                        annot["HGVSp"] = getConsistentHGVS(HGVS_by_tr[tr_base_acc]["HGVSp"], old_HGVSp)
                                        if old_HGVSp != annot["HGVSp"]:
                                            is_fixed_HGVSp = True
                                new_annot.append(annot)
                        record.info[args.annotations_field] = new_annot
                        # Trace results
                        if is_fixed_HGVSg:
                            nb_records["fixed_HGVSg"] += 1
                        if is_fixed_HGVSc:
                            nb_records["fixed_HGVSc"] += 1
                        if is_fixed_HGVSp:
                            nb_records["fixed_HGVSp"] += 1
                        if contains_colloc_annot:
                            nb_records["contains_colloc_annot"] += 1
                FH_out.write(record)
    # Log
    if nb_records["contains_colloc_annot"] != 0:
        log.warning("{}/{} variants contain collocated annotation removed during the process.".format(nb_records["contains_colloc_annot"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on HGVSg.".format(nb_records["fixed_HGVSg"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSc.".format(nb_records["fixed_HGVSc"], nb_records["analysed"]))
    log.info("{}/{} variants have been fixed on one of their HGVSp.".format(nb_records["fixed_HGVSp"], nb_records["analysed"]))
    log.info("End of job")
