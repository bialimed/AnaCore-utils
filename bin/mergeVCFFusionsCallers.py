#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.2'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import numpy
import logging
import argparse
import itertools
from anacore.vcf import VCFIO, HeaderInfoAttr, HeaderFormatAttr, getAlleleRecord
from anacore.region import Region, RegionList, iterOverlappedByRegion


########################################################################
#
# FUNCTIONS
#
########################################################################
def getNewHeaderAttr(args):
    """
    Return renamed and new VCFHeader elements for the merged VCF.

    :param args: The script's parameters.
    :type args: NameSpace
    :return: VCFHeader elements (filter, info, format, samples).
    :rtype: dict
    """
    unchanged_info = {"MATEID", "RNA_FIRST", "SVTYPE"}
    final_filter = {}
    final_info = {
        "CIPOS": HeaderInfoAttr(
            "CIPOS", type="Integer", number="2", description="Confidence interval around POS"
        ),
        "IDSRC": HeaderInfoAttr(
            "IDSRC", type="String", number=".", description="ID of breakend by source"
        ),
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


def getBNDInterval(record):
    """
    Return the start and end for BND record (CIPOS is taken into account).

    :param record: The breakend record.
    :type record: anacore.vcf.VCFRecord
    :return: Start and end for BND record (CIPOS is taken into account).
    :rtype: (int, int)
    """
    start = record.pos
    end = record.pos
    if "CIPOS" in record.info:
        start += record.info["CIPOS"][0]
        end += record.info["CIPOS"][1]
        if "IMPRECISE" in record.info:
            start -= 1
    return start, end


def getCount(data, field):
    """
    Return the value contained in the field specified for the first alternative variant.

    :param data: Values contained in sample field of a VCFRecord.
    :type data: dict
    :param field: Field title for the count.
    :type field: str
    :return: The value contained in the field specified for the first alternative variant.
    :rtype: int
    """
    val = 0
    if field in data:
        if not isinstance(data[field], list):
            val = data[field]
        else:
            if len(data[field]) == 1:
                val = data[field][0]
            else:
                val = data[field][1]
    return val


def loadBNDByID(in_vcf):
    """
    Return breakend by ID from a VCF file.

    :param in_vcf: Path to the VCF containing BND coming from one fusion caller (format: VCF).
    :type in_vcf: str
    :return: Breakend by ID.
    :rtype: dict
    """
    bnd_by_id = {}
    with VCFIO(in_vcf) as reader:
        if "SR" in reader.info and reader.info["SR"].number == ".":
            raise Exception('The number attribute for SR must be "A" or "R" or "1".')
        if "PR" in reader.info and reader.info["PR"].number == ".":
            raise Exception('The number attribute for PR must be "A" or "R" or "1".')
        for record in reader:
            if record.info["SVTYPE"] == "BND":
                bnd_by_id[record.id] = record
    return bnd_by_id


def groupBNDByFusions(bnd_by_id, annotation_field):
    """
    Return by chromosome the region of the first breakend in each fucion. The annotation of regions contains the two breakends (tags: first and second).

    :param bnd_by_id: Breakend by ID coming from one fusion caller.
    :type bnd_by_id: dict
    :param annotation_field: Field used to store annotations.
    :type annotation_field: str
    :return: By chromosome the region of the first breakend in each fucion. The annotation of regions contains the two breakends (tags: first and second).
    :rtype: dict
    """
    caller_fusions = dict()
    processed_fusions = set()
    fusion_by_name = {}
    for id, record in bnd_by_id.items():
        for alt_idx, alt in enumerate(record.alt):
            alt_first_bnd = record
            first_new_id = alt_first_bnd.id
            if len(record.alt) > 1:
                first_new_id += "_" + str(alt_idx)  # Record must be splitted for each mate
                alt_first_bnd = getAlleleRecord(record, alt_idx)
                alt_first_bnd.info["MATEID"] = [record.info["MATEID"][alt_idx]]
            mate_id = alt_first_bnd.info["MATEID"][0]
            mate_record = bnd_by_id[mate_id]
            alt_second_bnd = mate_record
            second_new_id = alt_second_bnd.id
            if len(mate_record.alt) > 1:
                first_idx = mate_record.info["MATEID"].index(alt_first_bnd.id)
                second_new_id += "_" + first_idx  # Record must be splitted for each mate
                alt_second_bnd = getAlleleRecord(mate_record, first_idx)
                alt_second_bnd.info["MATEID"] = [mate_record.info["MATEID"][first_idx]]
            fusion_id = " @@ ".join(sorted([alt_first_bnd.id, alt_second_bnd.id]))
            alt_first_bnd.id = first_new_id
            alt_second_bnd.info["MATEID"] = [first_new_id]
            alt_second_bnd.id = second_new_id
            alt_first_bnd.info["MATEID"] = [second_new_id]
            if fusion_id not in processed_fusions:
                processed_fusions.add(fusion_id)
                if "RNA_FIRST" not in alt_first_bnd.info and "RNA_FIRST" not in alt_second_bnd.info:
                    raise Exception("Tag RNA_FIRST must be present in one of the breakend {} or {}.".format(alt_first_bnd.id, mate_id))
                if "RNA_FIRST" in alt_second_bnd.info:
                    aux = alt_first_bnd
                    alt_first_bnd = alt_second_bnd
                    alt_second_bnd = aux
                interval_first_bnd = getBNDInterval(alt_first_bnd)
                fusion_name = " @@ ".join(sorted([alt_first_bnd.getName(), alt_second_bnd.getName()]))
                if fusion_name not in fusion_by_name:
                    region_first_bnd = Region(
                        interval_first_bnd[0],
                        interval_first_bnd[1],
                        reference=alt_first_bnd.chrom,
                        annot={"first": alt_first_bnd, "second": alt_second_bnd}
                    )
                    if alt_first_bnd.chrom not in caller_fusions:
                        caller_fusions[alt_first_bnd.chrom] = RegionList()
                    caller_fusions[alt_first_bnd.chrom].append(region_first_bnd)
                    fusion_by_name[fusion_name] = region_first_bnd
                else:  # Caller contains several entries for the same pair of breakends (same fusion but several anotations)
                    fusion_by_name[fusion_name].annot["first"].info[annotation_field] += alt_first_bnd.info[annotation_field]
                    fusion_by_name[fusion_name].annot["second"].info[annotation_field] += alt_second_bnd.info[annotation_field]
    return caller_fusions


def renameFields(bnd_record, caller_prefix, shared_filters):
    """
    Rename fields with prefix of the variant caller.

    :param bnd_record: The breakend record.
    :type bnd_record: anacore.vcf.VCFRecord
    :param caller_prefix: Prefix added to field names.
    :type caller_prefix: str
    :param shared_filters: Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix.
    :type shared_filters: set
    """
    unchanged_info = {"MATEID", "RNA_FIRST", "SVTYPE"}
    # Rename filters
    if bnd_record.filter is not None:
        new_filter = []
        for tag in bnd_record.filter:
            if tag != "PASS":
                if tag in shared_filters:  # Rename filters not based on caller
                    new_filter.append(tag)
                else:
                    new_filter.append("{}_{}".format(caller_prefix, tag))
        bnd_record.filter = new_filter
    # Rename INFO
    new_info = {}
    for key, val in bnd_record.info.items():
        if key in unchanged_info:
            new_info[key] = val
        else:
            new_info["{}_{}".format(caller_prefix, key)] = val
    bnd_record.info = new_info
    # Backup quality
    if bnd_record.qual is not None:
        bnd_record.info["{}_VCQUAL".format(caller_prefix)] = bnd_record.qual
    # Rename FORMAT
    bnd_record.format = ["{}_{}".format(caller_prefix, curr_filter) for curr_filter in bnd_record.format]
    for spl_name, spl_info in bnd_record.samples.items():
        renamed_info = {}
        for key, val in spl_info.items():
            renamed_info["{}_{}".format(caller_prefix, key)] = val
        bnd_record.samples[spl_name] = renamed_info


def getPrevFusion(records_pair, overlapped, curr_caller):
    """
    Return fusion with same right breakpoint of evaluated fusion.

    :param records_pair: Pair of breakends for the evaluated fusion.
    :type records_pair: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
    :param overlapped: List of breakends pairs with first breakend overlapping the first breakend of the evaluated fusion.
    :type overlapped: list
    :param curr_caller: Name of the caller used to produce evaluated fusion.
    :type curr_caller: str
    :return: Fusion with same right breakpoint of evaluated fusion.
    :rtype: (anacore.vcf.VCFRecord, anacore.vcf.VCFRecord)
    """
    prev_records = None
    left_record, right_record = records_pair
    if len(overlapped) > 0:
        start_right_record, end_right_record = getBNDInterval(right_record)
        for overlap_eval in overlapped:
            start_right_eval, end_right_eval = getBNDInterval(overlap_eval.annot["second"])
            if not start_right_record > end_right_eval and not end_right_record < start_right_eval:
                if prev_records is not None:
                    raise Exception(
                        "{} from {} has an overlap ambiguity between: {} and {}.".format(
                            left_record.getName(),
                            curr_caller,
                            prev_records[0].getName(),
                            overlap_eval.annot["first"].getName()
                        )
                    )
                prev_records = (overlap_eval.annot["first"], overlap_eval.annot["second"])
                log.debug(
                    "Merge {} from {} with {} from {}.".format(
                        overlap_eval.annot["first"].getName(),
                        " and ".join(prev_records[0].info["SRC"]),
                        left_record.getName(),
                        curr_caller
                    )
                )
    return prev_records


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
    whole_fusions = {}  # fisrt bnd region by chromosome
    for idx_in, curr_in in enumerate(inputs_variants):
        curr_caller = calling_sources[idx_in]
        log.info("Process {}".format(curr_caller))
        # breakend by id
        bnd_by_id = loadBNDByID(curr_in)
        # Group by fusion
        curr_caller_fusions = groupBNDByFusions(bnd_by_id, annotation_field)
        # Merge to other callers
        new_fusions = []
        for chrom, query, overlapped in iterOverlappedByRegion(curr_caller_fusions, whole_fusions):
            records = (query.annot["first"], query.annot["second"])
            # Extract PR and SR
            support_by_spl = {}
            for spl, data in records[0].samples.items():
                support_by_spl[spl] = {
                    "PR": getCount(data, "PR"),
                    "SR": getCount(data, "SR")
                }
            # Rename fields
            for curr_record in records:
                renameFields(curr_record, "s{}".format(idx_in), shared_filters)
            # Add to storage
            prev_records = getPrevFusion(records, overlapped, curr_caller)  # Get identical fusion from previous callers
            if prev_records is None:  # Prepare new fusion
                new_fusions.append(query)
                for curr_record in records:
                    # Data source
                    curr_record.info["SRC"] = [curr_caller]
                    curr_record.info["IDSRC"] = [curr_record.id]
                    # CIPOS
                    if "s{}_CIPOS".format(idx_in) in curr_record.info:
                        curr_record.info["CIPOS"] = curr_record.info["s{}_CIPOS".format(idx_in)]
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
            else:  # Update previous fusion
                for prev_rec, curr_rec in zip(prev_records, records):
                    prev_rec.info["SRC"].append(curr_caller)
                    prev_rec.info["IDSRC"].append(curr_rec.id)
                    # FILTERS
                    if curr_rec.filter is not None:
                        if prev_rec.filter is None:
                            prev_rec.filter = curr_rec.filter
                        else:
                            prev_rec.filter = list(set(prev_rec.filter) or set(curr_rec.filter))
                    # FORMAT
                    prev_rec.format.extend(curr_rec.format)
                    # INFO
                    del(curr_rec.info["MATEID"])
                    prev_rec.info.update(curr_rec.info)
                    # SAMPLES
                    for spl_name, spl_data in prev_rec.samples.items():
                        spl_data.update(curr_rec.samples[spl_name])
                        spl_data["SRSRC"].append(support_by_spl[spl_name]["SR"])
                        spl_data["PRSRC"].append(support_by_spl[spl_name]["PR"])
        # Add new fusions in whole_fusions
        for curr in new_fusions:
            if curr.reference.name not in whole_fusions:
                whole_fusions[curr.reference.name] = RegionList()
            whole_fusions[curr.reference.name].append(curr)
        # Sort fusions by first breakend
        for chrom, fusions in whole_fusions.items():
            whole_fusions[chrom] = RegionList(sorted(fusions, key=lambda x: (x.start, x.end)))
    # Flatten fusions
    returned_fusions = []
    for chr, fusions in whole_fusions.items():
        for fusion_region in fusions:
            returned_fusions.append((
                fusion_region.annot["first"],
                fusion_region.annot["second"]
            ))
    return returned_fusions


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
            log.info("Differences between retained {} and others callers (without missing): median={}, upper_quartile={}, 90_persentile={} and max={} on {} variants".format(
                metric,
                numpy.percentile(diff_by_metric[metric], 50, interpolation='midpoint'),
                numpy.percentile(diff_by_metric[metric], 75, interpolation='midpoint'),
                numpy.percentile(diff_by_metric[metric], 90, interpolation='midpoint'),
                max(diff_by_metric[metric]),
                nb_var
            ))


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
    parser = argparse.ArgumentParser(description='Merge VCF coming from different fusions caller on same sample(s). It is strongly recommended to apply this script before annotation and filtering/tagging.')
    parser.add_argument('-a', '--annotation-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-c', '--calling-sources', required=True, nargs='+', help='Name of the source in same order of --inputs-variants.')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-s', '--shared-filters', nargs='*', default=[], help='Filters tags applying to the variant and independent of caller like filters on annotations. These filters are not renamed to add caller ID as suffix. [Default: %(default)s]')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--inputs-variants', required=True, nargs='+', help='Path to the variants files coming from different callers (format: VCF). The order determine the which SR and PR are retained: the first caller where it is found in this list.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to the merged variants file (format: VCF).')
    args = parser.parse_args()
    args.shared_filters = set(args.shared_filters)

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
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
