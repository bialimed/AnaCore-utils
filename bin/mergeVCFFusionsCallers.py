#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import numpy
import logging
import argparse
import itertools
from anacore.vcf import VCFIO, HeaderInfoAttr, HeaderFormatAttr, getAlleleRecord


########################################################################
#
# FUNCTIONS
#
########################################################################
def cmpUniqName(record):
    alt = None
    if record.alt[0].startswith("[") or record.alt[0].startswith("]"):
        alt = "{}N".format(record.alt[0][:-1])
    else:
        alt = "N{}".format(record.alt[0][1:])
    return "{}:{}={}".format(record.chrom, record.pos, alt)


def getNewHeaderAttr(args):
    """
    Return renamed and new VCFHeader elements for the merged VCF.

    :param args: The script's parameters.
    :type args: NameSpace
    :return: VCFHeader elements (filter, info, format, samples).
    :rtype: dict
    """
    unchanged_info = {"MATEID", "RNA_FIRST", "SVTYPE", args.annotation_field}
    final_filter = {}
    final_info = {
        "SRC": HeaderInfoAttr(
            "SRC", type="String", number=".", description="Fusions callers where the breakend is identified. Possible values: {}".format(
                {name: "s" + str(idx) for idx, name in enumerate(args.calling_sources)}
            )
        )
    }
    final_format = {
        "SR": HeaderFormatAttr("SR", type="Integer", number="1", description="Count of reads mapping on the fusion junction"),
        "PR": HeaderFormatAttr("PR", type="Integer", number="1", description="Count of pairs of reads supporting the fusion"),
        "SRSRC": HeaderFormatAttr("SRSRC", type="Integer", number=".", description="Count of reads mapping on the fusion junction by source"),
        "PRSRC": HeaderFormatAttr("PRSRC", type="Integer", number=".", description="Count of pairs of reads supporting the fusion by source")
    }
    final_samples = None
    for idx_in, curr_in in enumerate(args.inputs_variants):
        with VCFIO(curr_in) as FH_vcf:
            # Samples
            if final_samples is None:
                final_samples = FH_vcf.samples
            elif FH_vcf.samples != final_samples:
                raise Exception(
                    "The samples in VCF are not the same: {} in {} and {} in {}.".format(
                        final_samples,
                        args.inputs_variants[0],
                        FH_vcf.samples,
                        curr_in
                    )
                )
            # FILTER
            for tag, data in FH_vcf.filter.items():
                new_tag = tag
                if tag not in args.shared_filters:  # Rename filters not based on caller
                    new_tag = "s{}_{}".format(idx_in, tag)
                    data.id = new_tag
                    data.source = args.calling_sources[idx_in]
                final_filter[new_tag] = data
            # INFO
            for tag, data in FH_vcf.info.items():
                if tag in unchanged_info:
                    if tag not in final_info or len(final_info[tag].description) < len(data.description):  # Manage merge between callers with 0 variants (and 0 annotations) and callers with variants
                        final_info[tag] = data
                else:
                    new_tag = "s{}_{}".format(idx_in, tag)
                    data.id = new_tag
                    data.source = args.calling_sources[idx_in]
                    final_info[new_tag] = data
            qual_tag = "s{}_VCQUAL".format(idx_in)
            final_info[qual_tag] = HeaderInfoAttr(qual_tag, type="Float", number="1", description="The variant quality", source=args.calling_sources[idx_in])
            # FORMAT
            for tag, data in FH_vcf.format.items():
                new_tag = "s{}_{}".format(idx_in, tag)
                data.id = new_tag
                data.source = args.calling_sources[idx_in]
                final_format[new_tag] = data
    return {
        "filter": final_filter,
        "info": final_info,
        "format": final_format,
        "samples": final_samples
    }


def getMergedRecords(inputs_variants, calling_sources, annotation_field, shared_filters):
    """
    Merge VCFRecords coming from several variant callers.

    :param inputs_variants: Pathes to the variants files.
    :type inputs_variants: list
    :param calling_sources: Names of the variants callers (in same order as inputs_variants).
    :type calling_sources: list
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :param shared_filters: Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix.
    :type shared_filters: set
    :return: Merged VCF records.
    :rtype: list
    """
    unchanged_info = {"MATEID", "RNA_FIRST", "SVTYPE", annotation_field}
    whole_by_name = {}
    for idx_in, curr_in in enumerate(inputs_variants):
        curr_caller = calling_sources[idx_in]
        log.info("Process {}".format(curr_caller))
        # breakend by id
        bnd_by_id = {}
        with VCFIO(curr_in) as reader:
            for record in reader:
                bnd_by_id[record.id] = record
        # Group by fusion
        fusion_by_name = {}
        for id, record in bnd_by_id.items():
            for alt_idx, alt in enumerate(record.alt):
                alt_first_bnd = record
                if len(record.alt) > 1:
                    alt_first_bnd = getAlleleRecord(record, alt_idx)
                    alt_first_bnd.info["MATEID"] = [record.info["MATEID"][alt_idx]]
                mate_id = alt_first_bnd.info["MATEID"][0]
                mate_record = bnd_by_id[mate_id]
                alt_second_bnd = mate_record
                if len(mate_record.alt) > 1:
                    first_idx = mate_record.info["MATEID"].index(alt_first_bnd.id)
                    alt_second_bnd = getAlleleRecord(mate_record, first_idx)
                    alt_second_bnd.info["MATEID"] = [mate_record.info["MATEID"][first_idx]]
                if "RNA_FIRST" not in alt_first_bnd.info and "RNA_FIRST" not in alt_second_bnd.info:
                    raise Exception("Tag RNA_FIRST must be present in one of the breakend {} or {}.".format(alt_first_bnd.id, mate_id))
                if "RNA_FIRST" in alt_second_bnd.info:
                    aux = alt_first_bnd
                    alt_first_bnd = alt_second_bnd
                    alt_second_bnd = aux
                fusion_by_name[cmpUniqName(alt_first_bnd) + "@" + cmpUniqName(alt_second_bnd)] = (alt_first_bnd, alt_second_bnd)
        # Merge to other callers
        for fusion_name, records in fusion_by_name.items():
            # Extract PR and SR
            support_by_spl = {}
            for spl, data in records[0].samples.items():
                support_by_spl[spl] = {
                    "PR": data["PR"] if "PR" in data else 0,
                    "SR": data["SR"] if "SR" in data else 0
                }
            # Rename fields
            for curr_record in records:
                # Rename filters
                if curr_record.filter is not None:
                    new_filter = []
                    for tag in curr_record.filter:
                        if tag != "PASS":
                            if tag in shared_filters:  # Rename filters not based on caller
                                new_filter.append(tag)
                            else:
                                new_filter.append("s{}_{}".format(idx_in, tag))
                    curr_record.filter = new_filter
                # Rename INFO
                new_info = {}
                for key, val in curr_record.info.items():
                    if key in unchanged_info:
                        new_info[key] = val
                    else:
                        new_info["s{}_{}".format(idx_in, key)] = val
                curr_record.info = new_info
                # Backup quality
                if curr_record.qual is not None:
                    curr_record.info["s{}_VCQUAL".format(idx_in)] = curr_record.qual
                # Rename FORMAT
                curr_record.format = ["s{}_{}".format(idx_in, curr_filter) for curr_filter in curr_record.format]
                for spl_name, spl_info in curr_record.samples.items():
                    renamed_info = {}
                    for key, val in spl_info.items():
                        renamed_info["s{}_{}".format(idx_in, key)] = val
                    curr_record.samples[spl_name] = renamed_info
            # Add to storage
            left_record, right_record = records
            if fusion_name not in whole_by_name:
                whole_by_name[fusion_name] = records
                for curr_record in records:
                    # Data source
                    curr_record.info["SRC"] = [curr_caller]
                    # Quality
                    if idx_in != 0:
                        curr_record.qual = None  # For consistency, the quality of the variant comes only from the first caller of the variant
                    # SR and PR by sample (from the first caller finding the variant: callers are in user order)
                    curr_record.format.insert(0, "SRSRC")
                    curr_record.format.insert(0, "PRSRC")
                    curr_record.format.insert(0, "SR")
                    curr_record.format.insert(0, "PR")
                    for spl_name, spl_data in curr_record.samples.items():
                        spl_data["SR"] = support_by_spl[spl_name]["SR"]
                        spl_data["PR"] = support_by_spl[spl_name]["PR"]
                        spl_data["SRSRC"] = [support_by_spl[spl_name]["SR"]]
                        spl_data["PRSRC"] = [support_by_spl[spl_name]["PR"]]
            else:
                prev_records = whole_by_name[fusion_name]
                for prev_rec, curr_rec in zip(prev_records, records):
                    prev_rec.info["SRC"].append(curr_caller)
                    # IDs
                    if curr_rec.id is not None:
                        prev_ids = prev_rec.id.split(";")
                        prev_ids.extend(curr_rec.id.split(";"))
                        prev_ids = sorted(list(set(prev_ids)))
                        prev_rec.id = ";".join(prev_ids)
                    # FILTERS
                    if curr_rec.filter is not None:
                        if prev_rec.filter is None:
                            prev_rec.filter = curr_rec.filter
                        else:
                            prev_rec.filter = list(set(prev_rec.filter) or set(curr_rec.filter))
                    # FORMAT
                    prev_rec.format.extend(curr_rec.format)
                    # INFO
                    prev_rec.info.update(curr_rec.info)
                    for spl_name, spl_data in prev_rec.samples.items():
                        spl_data.update(curr_rec.samples[spl_name])
                        spl_data["SRSRC"].append(support_by_spl[spl_name]["SR"])
                        spl_data["PRSRC"].append(support_by_spl[spl_name]["PR"])
    return whole_by_name.values()


def logSupportVariance(fusions, log):
    """
    Display in log the variance on support counts (SR and PR) between callers.

    :param breakends: Merged VCF records.
    :type breakends: list
    :param log: Logger object.
    :type log: logging.Logger
    """
    diff_by_metric = {"SR": [], "PR": []}
    nb_var = 0
    for records in fusions:
        if len(records[0].info["SRC"]) > 1:
            nb_var += 1
            for spl_name, spl_data in records[0].samples.items():
                for metric in ["SR", "PR"]:
                    retained = spl_data[metric]
                    max_diff = 0
                    for curr in spl_data[metric + "SRC"][1:]:
                        max_diff = max(max_diff, abs(retained - curr))
                    diff_by_metric[metric].append(max_diff)
    # Log
    for metric in ["SR", "PR"]:
        if nb_var == 0:
            log.info("Differences between retained {} and others callers (without missing): 0 common variants".format(metric))
        else:
            log.info("Differences between retained {} and others callers (without missing): median={:.1%}, upper_quartile={:.1%}, 90_persentile={:.1%} and max={:.1%} on {} variants".format(
                metric,
                numpy.percentile(diff_by_metric[metric], 50, interpolation='midpoint'),
                numpy.percentile(diff_by_metric[metric], 75, interpolation='midpoint'),
                numpy.percentile(diff_by_metric[metric], 90, interpolation='midpoint'),
                max(diff_by_metric[metric]),
                nb_var
            ))


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Merge VCF coming from different fusions caller on same sample(s). It is strongly recommended to apply this script before annotation and filtering/tagging.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-s', '--shared-filters', nargs='*', default=[], help='Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-sources', required=True, nargs='+', help='Name of the source in same order of --inputs-variants.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', required=True, nargs='+', help='Path to the variants files coming from different callers (format: VCF). The order determine the which SR and PR are retained: the first caller where it is found in this list.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()
    args.shared_filters = set(args.shared_filters)

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get merged records
    fusions = getMergedRecords(args.inputs_variants, args.calling_sources, args.annotation_field, args.shared_filters)

    # Log differences in SR and PR
    logSupportVariance(fusions, log)

    # Write
    breakends = list(itertools.chain.from_iterable(fusions))
    with VCFIO(args.output_variants, "w") as writer:
        # Header
        new_header = getNewHeaderAttr(args)
        writer.samples = new_header["samples"]
        writer.info = new_header["info"]
        writer.format = new_header["format"]
        writer.filter = new_header["filter"]
        writer.writeHeader()
        # Records
        for record in sorted(breakends, key=lambda record: (record.chrom, record.refStart(), record.refEnd())):
            if record.filter is not None and len(record.filter) == 0:
                record.filter = ["PASS"]
            writer.write(record)

    log.info("End of job")
