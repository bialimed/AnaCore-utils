#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '2.3.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import pysam
import logging
import argparse
from copy import deepcopy
from statistics import mean
from collections import deque
from anacore.sequenceIO import IdxFastaIO
from anacore.vcf import VCFIO, VCFRecord, HeaderInfoAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def getAlnCmp(read, ref_seq):
    """
    Return the dense representation of the alignment between reference sequence and read.

    :param read: The alignment.
    :type read: pysam.AlignedSegment
    :param ref_seq: The reference sequence included in the alignement.
    :type ref_seq: str
    :return: Dense representation of the alignment. First is reference alignment and second is read alignment.
    :rtype: (list, list)
    """
    ref_aln = []
    read_aln = []
    ref_seq = [elt for elt in ref_seq]  # reference sequence on alignment (read.get_reference_sequence() can be incorrect)
    read_seq = [elt for elt in read.query_alignment_sequence]  # query sequence without clipped
    for operation_id, operation_lg in read.cigartuples:
        if operation_id != 5 and operation_id != 4:  # Is not clipped
            if operation_id in {0, 7, 8}:  # Match or mismatch
                read_aln.extend(read_seq[:operation_lg])
                read_seq = read_seq[operation_lg:]
                ref_aln.extend(ref_seq[:operation_lg])
                ref_seq = ref_seq[operation_lg:]
            elif operation_id == 1:  # Insertion
                if len(read_aln) != 0:  # Insertion does not start alignment
                    read_aln[-1] += "".join(read_seq[:operation_lg])
                    # Insertion starting alignment should be clipping.
                    # In this situation reference nucleotids cannot be found and
                    # reference start is after the insertion.
                read_seq = read_seq[operation_lg:]
            elif operation_id == 2:  # Deletion
                read_aln.extend(["" for operation_pos in range(operation_lg)])
                ref_aln.extend(ref_seq[:operation_lg])
                ref_seq = ref_seq[operation_lg:]
            elif operation_id == 3:  # Refskip ########################### Are the introns in the variant ?
                # next nt in ref but not in current read
                ref_seq = ref_seq[operation_lg:]
                ref_aln.extend(["/" for operation_pos in range(operation_lg)])
                read_aln.extend(["/" for operation_pos in range(operation_lg)])
            elif operation_id == 9:  # Back (www.seqanswers.com/forums/showthread.php?t=34440)
                raise Exception("Parsing error on read {}. The management for the CIGAR operator B is not implemented.".format(read.query_name))
            # elif operation_id == 6:  # Padding
            #     pass
    return ref_aln, read_aln


def setRefPos(variant, seq_handler, padding=200):
    """
    Add start and end attributes in VCFRecord. For insertions the start is defined on the first position before the insertion and the end on the last position affected by the insertion.

    :param variant: The variant to update.
    :type variant: anacore.vcf.VCFRecord
    """
    if variant.ref == VCFRecord.getEmptyAlleleMarker() or variant.alt[0] == VCFRecord.getEmptyAlleleMarker():  # Normalized indel
        # Most upstream
        variant.upstream_start, variant.upstream_end = getStartEnd(variant)
        # Most downstream
        sub_region = seq_handler.getSub(variant.chrom, variant.pos - 2, variant.pos + len(variant.ref) + padding)
        chrom_pos = variant.pos
        variant.pos = 3  # Switch position from chromosome to position from subregion
        downstream_var = variant.getMostDownstream(sub_region)
        variant.pos = chrom_pos + variant.pos - 3  # Switch position from subregion to position from chromosome
        downstream_var.pos = variant.pos
        variant.downstream_start, variant.downstream_end = getStartEnd(downstream_var)
    else:
        variant.upstream_start, variant.upstream_end = getStartEnd(variant)
        variant.downstream_start = variant.upstream_start
        variant.downstream_end = variant.upstream_end


def getStartEnd(variant):
    """
    Return reference start and end for a variant. For insertions the start is defined on the first position before the insertion and the end on the last position affected by the insertion.

    :param variant: The variant.
    :type variant: anacore.vcf.VCFRecord
    :return: start and end.
    :rtype: (int, int)
    """
    start = int(variant.refStart())
    if variant.isInsertion():
        if start == int(variant.refStart() + 0.5):
            start -= 1
    end = int(variant.refEnd())
    return start, end


def getReadRefAlt(ref_aln, read_aln, ref_start, target_is_ins, target_start, target_end):
    """
    Return reference sequence and read sequence for the selected target.

    :param ref_aln: Reference sequence in alignment.
    :type ref_aln: str
    :param read_aln: Read sequence in alignment.
    :type read_aln: str
    :param ref_start: Start position for the alignment on the reference (1-based).
    :type ref_start: int
    :param target_is_ins: The anlyzed target correspond to an exon.
    :type target_is_ins: bool
    :param target_start: Start reference position for the target (1-based).
    :type target_start: int
    :param target_end: End reference position for the target (1-based).
    :type target_end: int
    :return: Reference and alternative sequence for the read.
    :rtype: (str, str)
    """
    alt = read_aln[target_start - ref_start:target_end - ref_start + 1]
    ref = ref_aln[target_start - ref_start:target_end - ref_start + 1]
    if target_is_ins:
        while len(ref[0]) > 0 and len(alt[0]) > 0 and alt[0][0] == ref[0][0]:
            alt[0] = alt[0][1:]
            ref[0] = ref[0][1:]
    return ref, alt


def getIncludingReadsDNA(FH_aln, FH_seq, chrom_id, target_start, target_end):
    """
    Return read ID of reads including the target.

    :param FH_aln: The file handle to the alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param FH_seq: File handle to the refersence sequence file. Unused in this function but add to keep the same interface with getIncludingReadsRNA.
    :type FH_seq: IdxFastaIO
    :param chrom_id: Chromosome ID.
    :type chrom_id: str
    :param target_start: Start position for target.
    :type target_start: int
    :param target_end: End position for target.
    :type target_end: int
    :return: Reads IDs of reads including the target.
    :rtype: set
    """
    including_reads = set()
    for read in FH_aln.fetch(chrom_id, target_start - 1, target_end):
        if not read.is_duplicate and not read.is_secondary:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                includes = (ref_start <= target_start and ref_end >= target_end)
                if includes:
                    including_reads.add(read.query_name)
    return including_reads


def getIncludingReadsRNA(FH_aln, FH_seq, chrom_id, target_start, target_end):
    """
    Return read ID of reads including the target.

    :param FH_aln: File handle to the alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param FH_seq: File handle to the refersence sequence file.
    :type FH_seq: IdxFastaIO
    :param chrom_id: Chromosome ID.
    :type chrom_id: str
    :param target_start: Start position for target.
    :type target_start: int
    :param target_end: End position for target.
    :type target_end: int
    :return: Reads IDs of reads including the target.
    :rtype: set
    """
    including_reads = set()
    for read in FH_aln.fetch(chrom_id, target_start - 1, target_end):
        if not read.is_duplicate and not read.is_secondary:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                includes = (ref_start <= target_start and ref_end >= target_end)
                if includes:
                    ref_aln, read_aln = getAlnCmp(read, FH_seq.getSub(chrom_id, ref_start, ref_end))
                    target_ref_aln = ref_aln[target_start - ref_start:target_end - ref_start + 1]
                    if "/" not in set(target_ref_aln):
                        including_reads.add(read.query_name)
    return including_reads


def getSupportingReads(var, FH_seq, chrom_id, FH_aln, log):
    """
    Return read ID of reads supporting the altenative variant.

    :param var: The variant.
    :type var: anacore.vcf.VCFRecord updated with iniVariant() and isIns
    :param FH_seq: File handle to the refersence sequence file.
    :type FH_seq: IdxFastaIO
    :param chrom_id: Chromosome ID.
    :type chrom_id: str
    :param FH_aln: The file handle to the alignments file. The variants must have been defined from this alignments file.
    :type FH_aln: pysam.AlignmentFile
    :param log: The logger object.
    :type log: logging.Logger
    :return: The list of supporting reads IDs.
    :rtype: set
    """
    supporting_reads = set()
    is_insertion = var.isInsertion()
    for read in FH_aln.fetch(var.chrom, var.upstream_start - 1, var.downstream_end):
        if not read.is_duplicate and not read.is_secondary:
            reads_pos = read.get_reference_positions()
            if len(reads_pos) != 0:  # Skip alignment with problem
                ref_start = reads_pos[0] + 1  # 0-based to 1-based
                ref_end = reads_pos[-1] + 1  # 0-based to 1-based
                overlap_var = (ref_start <= var.upstream_start and ref_end >= var.downstream_end)
                if overlap_var:
                    ref_aln, read_aln = getAlnCmp(read, FH_seq.getSub(chrom_id, ref_start, ref_end))
                    var_alt = var.alt[0].upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                    var_ref = var.ref.upper().replace(VCFRecord.getEmptyAlleleMarker(), "")
                    # Test with upstream coordinates
                    ref, alt = getReadRefAlt(ref_aln, read_aln, ref_start, is_insertion, var.upstream_start, var.upstream_end)
                    if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # The alternative is present on most upstream coordinates
                        log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                        supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping
                    # Test with downstream coordinates
                    elif var.upstream_start != var.downstream_start:
                        ref, alt = getReadRefAlt(ref_aln, read_aln, ref_start, is_insertion, var.downstream_start, var.downstream_end)
                        if "".join(alt).upper() == var_alt and "".join(ref).upper() == var_ref:  # The alternative is present on most downstream coordinates
                            log.debug("{}\t{}/{}\t'{}'\t'{}'\t{}".format(read.query_name, var.ref, var.alt[0], "".join(ref), "".join(alt), read.cigarstring))
                            supporting_reads.add(read.query_name)  # Fragment is overlapping if at least one of his read is ovelapping
    return supporting_reads


def areColocated(first, second):
    """
    Return True if one of the two variants is included in the location of the other.

    :param first: The variant.
    :type first: anacore.vcf.VCFRecord updated with iniVariant()
    :param second: The variant.
    :type second: anacore.vcf.VCFRecord updated with iniVariant()
    :return: True if one of the two variants is included in the location of the other.
    :rtype: bool
    """
    first_up = set([elt for elt in range(first.upstream_start, first.upstream_end + 1)])
    second_up = set([elt for elt in range(second.upstream_start, second.upstream_end + 1)])
    first_down = set([elt for elt in range(first.downstream_start, first.downstream_end + 1)])
    second_down = set([elt for elt in range(second.downstream_start, second.downstream_end + 1)])
    included = False
    if len(first_up - second_up) == 0:
        included = True
    elif len(first_up - second_down) == 0:
        included = True
    elif len(first_down - second_up) == 0:
        included = True
    elif len(first_down - second_down) == 0:
        included = True
    elif len(second_up - first_up) == 0:
        included = True
    elif len(second_down - first_up) == 0:
        included = True
    elif len(second_up - first_down) == 0:
        included = True
    elif len(second_down - first_down) == 0:
        included = True
    return included


def initVariant(record, FH_seq):
    """
    Normalize and add attributes (supporting_reads to None, upstream_start, upstream_end, downstream_start, downstream_end).

    :param record: Record to update.
    :type record: anacore.vcf.VCFRecord
    :param FH_seq: File handle to the reference sequences file.
    :type FH_seq: anacore.sequenceIO.IdxFastaIO
    """
    record.normalizeSingleAllele()
    record.supporting_reads = None
    setRefPos(record, FH_seq)


def traceMerge(record, intersection_rate, intersection_count):
    """
    Trace merge info in record.

    :param record: Record to update.
    :type record: anacore.vcf.VCFRecord
    :param intersection_rate: Number of reads supporting the two variants against the number of reads supporting the two or only one of the two variants.
    :type intersection_rate: float
    :param intersection_count: Number of reads supporting the two variants.
    :type intersection_count: int
    """
    record.info["MCO_IR"] = []
    if "MCO_IR" in prev.info:
        for curr_IR in prev.info["MCO_IR"]:
            record.info["MCO_IR"].append(curr_IR)
    record.info["MCO_IR"].append(round(intersection_rate, 4))
    record.info["MCO_IC"] = []
    if "MCO_IC" in prev.info:
        for curr_IC in prev.info["MCO_IC"]:
            record.info["MCO_IC"].append(curr_IC)
    record.info["MCO_IC"].append(intersection_count)


def mergedRecord(vcf, first, first_std_name, second, second_std_name, FH_seq):
    """
    Return the VCFRecord corresponding to the merge of first and second.

    :param vcf: The file handle to VCF.
    :type vcf: anacore.vcf.VCFIO
    :param first: The upstream variant to merge.
    :type first: anacore.vcf.VCFRecord
    :param first_std_name: The initial name of the upstream variant to merge (before normalisation).
    :type first_std_name: str
    :param second: The downstream variant to merge.
    :type second: anacore.vcf.VCFRecord
    :param second_std_name: The initial name of the downstream variant to merge (before normalisation).
    :type second_std_name: str
    :param FH_seq: File handle to the refersence sequence file.
    :type FH_seq: IdxFastaIO
    :return: The variant corresponding to the merge of first and second.
    :rtype: anacore.vcf.VCFRecord
    :todo: Keep INFO and format on strand from FreeBayes, VarDict, ...
    """
    merged = VCFRecord(
        first.chrom,  # chrom
        first.pos,  # pos
        pFormat=first.format
    )
    # Ref and Alt
    first_end = int(round(first.refEnd() - 0.49, 0))
    second_start = int(round(second.refStart() + 0.49, 0))
    ref_add = ""
    if second_start - first_end > 0:
        ref_add = FH_seq.getSub(first.chrom, first_end + 1, second_start - 1)
    merged.ref = first.ref + ref_add + second.ref
    merged.ref = merged.ref.replace(VCFRecord.getEmptyAlleleMarker(), "")
    merged.alt = [first.alt[0] + ref_add + second.alt[0]]
    merged.alt[0] = merged.alt[0].replace(VCFRecord.getEmptyAlleleMarker(), "")
    # Filter
    first_filters = [] if first.filter is None else first.filter
    second_filters = [] if second.filter is None else second.filter
    merged.filter = list(set(first_filters + second_filters))
    if len(merged.filter) > 1 and "PASS" in merged.filter:
        merged.filter.remove("PASS")
    # Samples
    for spl in first.samples:
        merged.samples[spl] = {}
        if "DP" in first.format:
            merged.samples[spl]["DP"] = min(
                first.getDP(spl), second.getDP(spl)
            )
        if "AD" in first.format:
            if vcf.format["AD"].number == "1":  # Contains one alt allele
                merged.samples[spl]["AD"] = min(first.samples[spl]["AD"], second.samples[spl]["AD"])
            else:
                merged.samples[spl]["AD"] = [min(first_AD, second_AD) for first_AD, second_AD in zip(first.samples[spl]["AD"], second.samples[spl]["AD"])]
        if "AF" in first.format:
            if vcf.format["AF"].number == "1":  # Contains one alt allele
                merged.samples[spl]["AF"] = min(first.samples[spl]["AF"], second.samples[spl]["AF"])
            else:
                merged.samples[spl]["AF"] = [min(first_AF, second_AF) for first_AF, second_AF in zip(first.samples[spl]["AF"], second.samples[spl]["AF"])]
    # INFO metrics
    if "AD" in first.info:
        if vcf.info["AD"].number == "1":  # Contains one alt allele
            merged.info["AD"] = merged.getPopAltAD()[0]
        elif vcf.info["AD"].number == "R":  # Contains ref and alt alleles
            merged.info["AD"] = [merged.getPopRefAD()] + merged.getPopAltAD()
        else:  # Contains only alt alleles
            merged.info["AD"] = merged.getPopAltAD()
    if "DP" in first.info:
        merged.info["DP"] = merged.getPopDP()
    if "AF" in first.info:
        if vcf.info["AF"].number == "1":  # Contains one alt allele
            merged.info["AF"] = merged.getPopAltAF()[0]
        elif vcf.info["AF"].number == "R":  # Contains ref and alt alleles
            merged.info["AF"] = [merged.getPopRefAF()] + merged.getPopAltAF()
        else:  # Contains only alt alleles
            merged.info["AF"] = merged.getPopAltAF()
    # INFO Parents
    merged.info["MCO_VAR"] = []
    if "MCO_VAR" in first.info:
        for parent in first.info["MCO_VAR"]:
            merged.info["MCO_VAR"].append(parent)
    else:
        merged.info["MCO_VAR"].append(first_std_name)
    if "MCO_VAR" in second.info:
        for parent in second.info["MCO_VAR"]:
            merged.info["MCO_VAR"].append(parent)
    else:
        merged.info["MCO_VAR"].append(second_std_name)
    # Quality
    merged.info["MCO_QUAL"] = []
    if "MCO_QUAL" in first.info:
        for qual in first.info["MCO_QUAL"]:
            merged.info["MCO_QUAL"].append(qual)
    else:
        merged.info["MCO_QUAL"].append(first.qual)
    if "MCO_QUAL" in second.info:
        for qual in second.info["MCO_QUAL"]:
            merged.info["MCO_QUAL"].append(qual)
    else:
        merged.info["MCO_QUAL"].append(second.qual)
    if None not in merged.info["MCO_QUAL"]:
        merged.qual = mean(merged.info["MCO_QUAL"])
    # Return
    return merged


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
    """
    :todo: count fragment overlapping the two and not the reads (increase length)
    """
    # Manage parameters
    parser = argparse.ArgumentParser(description='Groups variants occuring in same reads.')
    parser.add_argument('-d', '--max-distance', default=10, type=int, help='Maximum distance between two merged variants. [Default: %(default)s]')
    parser.add_argument('-f', '--AF-diff-rate', default=0.2, type=float, help='Maximum difference rate between AF of two merged variants. [Default: %(default)s]')
    parser.add_argument('-l', '--logging-level', default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], action=LoggerAction, help='The logger level. [Default: %(default)s]')
    parser.add_argument('-n', '--intersection-count', default=3, type=int, help='Minimum number of reads containing co-occurancy. [Default: %(default)s]')
    parser.add_argument('-p', '--spliced-aln', action='store_true', help='Use this option to manage spliced alignment.')
    parser.add_argument('-r', '--intersection-rate', default=0.9, type=float, help='Minimum ratio of co-occurancy (nb_reads_containing_the_two_variants / nb_reads_overlapping_the_two_variants_but_containing_only_one). [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-aln', required=True, help='Path to the alignment file (format: BAM).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the variants file (format: VCF). Variants must be ordered by position and should be move to upstream.')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to the reference sequences file (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the variant file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(args.logging_level)
    log.info("Command: " + " ".join(sys.argv))

    # Merge variants
    getIncludingReads = getIncludingReadsRNA if args.spliced_aln else getIncludingReadsDNA
    with IdxFastaIO(args.input_sequences) as FH_seq:
        with VCFIO(args.output_variants, "w") as FH_out:
            with pysam.AlignmentFile(args.input_aln, "rb") as FH_aln:
                with VCFIO(args.input_variants) as FH_vcf:
                    # Header
                    FH_out.copyHeader(FH_vcf)
                    FH_out.info["MCO_VAR"] = HeaderInfoAttr("MCO_VAR", "Name of the variants merged because their occur on same reads.", type="String", number=".")
                    FH_out.info["MCO_QUAL"] = HeaderInfoAttr("MCO_QUAL", "Qualities of the variants merged because their occur on same reads.", type="String", number=".")
                    FH_out.info["MCO_IR"] = HeaderInfoAttr("MCO_IR", "Co-occurancy rate between pairs of variants.", type="String", number=".")
                    FH_out.info["MCO_IC"] = HeaderInfoAttr("MCO_IC", "Co-occurancy count between pairs of variants.", type="String", number=".")
                    FH_out.writeHeader()
                    # Records
                    prev_chrom = None
                    chrom_var = deque()
                    prev_list = list()
                    for curr in FH_vcf:
                        std_curr = deepcopy(curr)
                        initVariant(curr, FH_seq)
                        merged_idx = set()
                        removed_idx = None
                        if prev_chrom is not None and prev_chrom != curr.chrom:  # The chromosome has change between previous and current variant
                            for prev, std_prev in prev_list:
                                chrom_var.append(std_prev)
                            for rec in sorted(chrom_var, key=lambda x: (x.pos, x.refEnd(), x.alt[0])):
                                FH_out.write(rec)
                            chrom_var = deque()
                            prev_chrom = curr.chrom
                        else:  # The chromosome is the same between previous and current variant
                            for idx, (prev, std_prev) in enumerate(prev_list[::-1]):
                                if removed_idx is None:  # Distance between current variant and previous most upstream variants is ok
                                    variants_distance = max(
                                        max(0, curr.upstream_start - prev.downstream_end),  # prev is before
                                        max(0, prev.upstream_start - curr.downstream_end)  # prev is after
                                    )
                                    if variants_distance > args.max_distance:  # The two records are too far
                                        removed_idx = len(prev_list) - 1 - idx
                                    elif areColocated(curr, prev):  # The two records are colocated
                                        log.debug("Skip colocated variants {} and {}.".format(prev.getName(), curr.getName()))
                                    else:  # The two records are close together
                                        prev_AF = prev.getPopAltAF()[0]
                                        curr_AF = curr.getPopAltAF()[0]
                                        AF_diff = 1 - (min(prev_AF, curr_AF) / max(prev_AF, curr_AF))
                                        log.debug("Allelels frequencies for {} and {}: {:.1%} and {:.1%} (diff rate: {:.2}).".format(prev.getName(), curr.getName(), prev_AF, curr_AF, AF_diff))
                                        if AF_diff <= args.AF_diff_rate:  # The two records have similar frequencies
                                            # Set supporting reads
                                            if prev.supporting_reads is None:
                                                prev.supporting_reads = getSupportingReads(prev, FH_seq, curr.chrom, FH_aln, log)
                                            if curr.supporting_reads is None:
                                                curr.supporting_reads = getSupportingReads(curr, FH_seq, curr.chrom, FH_aln, log)
                                            shared_reads = getIncludingReads(FH_aln, FH_seq, curr.chrom, min(prev.upstream_start, curr.upstream_start), max(prev.downstream_end, curr.downstream_end))
                                            # Check co-occurence
                                            if len(shared_reads) == 0:
                                                log.warning("No read overlap the two evaluated variants: {} and {}. In this condition the merge cannot be evaluated.".format(prev.getName(), curr.getName()))
                                            else:
                                                prev_support_shared = (prev.supporting_reads & shared_reads)
                                                curr_support_shared = (curr.supporting_reads & shared_reads)
                                                intersection_count = len(prev_support_shared & curr_support_shared)
                                                analysed_count = len(prev_support_shared | curr_support_shared)
                                                intersection_rate = 0.0 if analysed_count == 0 else intersection_count / analysed_count
                                                log.debug("{} and {} intersection rate: {:.5} ; number: {}.".format(prev.getName(), curr.getName(), intersection_rate, intersection_count))
                                                if intersection_rate >= args.intersection_rate and intersection_count >= args.intersection_count:
                                                    # Merge variants
                                                    first = prev
                                                    first_std_name = std_prev.getName()
                                                    second = curr
                                                    second_std_name = std_curr.getName()
                                                    if first.upstream_start > second.upstream_start:
                                                        first = curr
                                                        first_std_name = std_curr.getName()
                                                        second = prev
                                                        second_std_name = std_prev.getName()
                                                    merged = mergedRecord(FH_vcf, first, first_std_name, second, second_std_name, FH_seq)
                                                    traceMerge(merged, intersection_rate, intersection_count)
                                                    log.info("Merge {} and {} in {} (intersection: {:.2f} on {}]).".format(
                                                        prev.getName(), curr.getName(), merged.getName(), intersection_rate, analysed_count
                                                    ))
                                                    # Prepare merged to become prev
                                                    merged.fastStandardize(FH_seq, 200)
                                                    std_curr = deepcopy(merged)
                                                    curr = merged
                                                    initVariant(curr, FH_seq)
                                                    merged_idx.add(len(prev_list) - 1 - idx)
                            # Store the far records and remove them from the previous ones
                            if removed_idx is None:  # All the variants was close
                                # Remove individual version of merged variants
                                for idx in sorted(merged_idx, reverse=True):
                                    del(prev_list[idx])
                            else:  # Some variants was too far
                                # Remove individual version of merged variants
                                for idx in sorted(merged_idx, reverse=True):
                                    if idx > removed_idx:
                                        del(prev_list[idx])
                                # Push too far vairants in chrom_var
                                for idx in range(removed_idx + 1):
                                    record, std_record = prev_list.pop()
                                    if idx not in merged_idx:
                                        chrom_var.append(std_record)
                        prev_list.append((curr, std_curr))
                        prev_list = sorted(prev_list, key=lambda var: var[0].downstream_end)
                    # Last chromosome
                    for prev, std_prev in prev_list:
                        chrom_var.append(std_prev)
                    for rec in sorted(chrom_var, key=lambda x: (x.pos, x.refEnd(), x.alt[0])):
                        FH_out.write(rec)
    log.info("End of job")
