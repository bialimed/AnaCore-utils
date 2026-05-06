#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'

import os
import sys
import logging
import argparse
from anacore.vcf import HeaderInfoAttr
from anacore.fusion import BreakendVCFIO
from anacore.sequenceIO import IdxFastaIO
import re

# Opposite shard same stream
#          ---->
#         +++++++|           chr1 7 . G G[chr2:8[
#    chr1 NNNTAGG|NNNNNNN
#
#    chr2 NNNTAGG|NNNNNNN
#                |+++++++    chr2 8 . N ]chr1:7]N
#                  ---->
#
#    res: NNNTAGGG|NNNNNN or NNNTAGG|GNNNNNN or ... or NNN|TAGGGNNNNNN
#
# Opposite shard opposite stream
#          ---->
#         +++++++|            chr1 7 . G G]chr2:8]
#    chr1 NNNTAGG|NNNNNNN
#
#    chr2 NNNATCC|NNNNNNN
#         NNNTAGG|NNNNNNN
#                |+++++++     chr2 8 . N [chr1:7[N
#                  <----
#
#    res: NNNTAGGG|NNNNNN or NNNTAGG|GNNNNNN or ... or NNN|TAGGGNNNNNN
#
# Same shard opposite stream
#          ---->
#         +++++++|            chr1 7 . N N]chr2:7]
#    chr1 NNNNNNN|TAGGNNN
#
#    chr2 NNNCCTA|NNNNNNN
#         NNNGGAT|NNNNNNN
#         +++++++|            chr2 7 . A A]chr1:7]
#          <----
#
#    res: NNNNNNN|TAGGNNN or NNNNNNNT|AGGNNN or ... or NNNNNNNTAGG|NNN
#
# Same shard same stream
#          ---->
#         +++++++|            chr1 7 . N N[chr2:7[
#    chr1 NNNNNNN|TAGGNNN
#
#    chr2 NNNGGAT|NNNNNNN
#         +++++++|            chr2 7 . T T[chr1:7[
#          ---->
#
#    res: NNNNNNN|TAGGNNN or NNNNNNNT|AGGNNN or ... or NNNNNNNTAGG|NNN


########################################################################
#
# CUSTOMS BEFORE anacore 3.3.0
#
########################################################################
def decodedAlt(alt_str):
    """
    Return alternative information as dict from alternative string.

    :param alt_str: One alternative from VCF alternative field (example: A]chr1:45874121]).
    :type alt_str: str
    :return: Alternative information as dict from alternative string.
    :rtype: dict
    """
    start = alt_str[0]
    if start == "[" or start == "]":
        regexp = r"(.)(.+):(.+)[\[\]](.+)"
        bracket, chrom, pos, nt = re.fullmatch(regexp, alt_str).groups()
        reverse = bracket == "["
        ref_shard = "down"
    else:
        regexp = r"(.+)([\[\]])(.+):(.+)."
        nt, bracket, chrom, pos = re.fullmatch(regexp, alt_str).groups()
        reverse = bracket == "]"
        ref_shard = "up"
    return {"ref_shard": ref_shard, "reverse": reverse, "chrom": chrom, "pos": int(pos), "nt": nt}


def encodedAlt(alt):
    """
    Return allele string in VCF alternative field format (example: A]chr1:45874121]) from alternative allele dict.

    :param alt: One alternative in dict format (example: {"ref_shard": "up", "reverse": True, "chrom": "chr1", "pos": 45874121, "nt": 1}).
    :type alt: dict
    :return: Allele string in VCF alternative field format (example: A]chr1:45874121]) from alternative allele dict.
    :rtype: str
    """
    if alt["ref_shard"] == "up":
        bracket = "]" if alt["reverse"] else "["
        res = "{}{}{}:{}{}".format(alt["nt"], bracket, alt["chrom"], alt["pos"], bracket)
    else:
        bracket = "[" if alt["reverse"] else "]"
        res = "{}{}:{}{}{}".format(bracket, alt["chrom"], alt["pos"], bracket, alt["nt"])
    return res


########################################################################
#
# FUNCTIONS
#
########################################################################
dna_complement = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N',
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a', 'n': 'n'
}


def getComplement(seq):
    """
    Return complementary sequence.

    :param seq: The sequence.
    :type seq: str
    :return: The complementary sequence.
    :rtype: str
    """
    return "".join([dna_complement[nt] for nt in seq])


def replaceUnknownNt(record, seq_handler):
    """
    Replace undetermined ref and alt by the convenient nucleotid.

    :param record: The breakend record.
    :type record: anacore.vcf.VCFRecord
    :param seq_handler: Indexed reader for the reference genome used in fusion calling.
    :type seq_handler: anacore.sequenceIO.IdxFastaIO
    """
    if record.ref == "N":
        record.ref = seq_handler.getSub(record.chrom, record.pos, record.pos)
    for curr_idx, curr_alt in enumerate(record.alt):
        if curr_alt.startswith("N"):  # N[1:1008754[ or NATTCAT[1:1008754[ or N]1:1008754] or NATTCAT]1:1008754]
            record.alt[curr_idx] = record.ref + curr_alt[1:]
        elif curr_alt.endswith("N"):  # [1:1008754[N or [1:1008754[ATTCATN or ]1:1008754]N or ]1:1008754]ATTCATN
            record.alt[curr_idx] = curr_alt[:-1] + record.ref


def fastStandardize(first, second, seq_handler, padding=50):
    """
    Each breakend of the pair is placed at the left most position, and the uncertainty is represented with the CIPOS tag. The ALT string is then constructed assuming this choice.

    :param first: The breakend of the first shard in fusion (donor).
    :type first: anacore.vcf.VCFRecord
    :param second: The breakend of the second shard in fusion (acceptor).
    :type second: anacore.vcf.VCFRecord
    :param seq_handler: Indexed reader for the reference genome used in fusion calling.
    :type seq_handler: anacore.sequenceIO.IdxFastaIO
    :param padding: Number of nucleotids to inspect before and after the breakends: upstream and downstream movements are limited to this number of nucleotids.
    :type padding: int
    """
    first_alt = decodedAlt(first.alt[0])
    second_alt = decodedAlt(second.alt[0])
    before_first = seq_handler.getSub(first.chrom, max(first.pos - padding, 1), first.pos)
    before_second = seq_handler.getSub(second.chrom, max(first_alt["pos"] - padding, 1), first_alt["pos"])
    after_first = seq_handler.getSub(first.chrom, first.pos, first.pos + padding)
    after_second = seq_handler.getSub(second.chrom, first_alt["pos"], first_alt["pos"] + padding)
    if first_alt["reverse"]:
        before_second = getComplement(before_second)
        after_second = getComplement(after_second)
    cipos_start = 0
    cipos_end = 0
    if first_alt["ref_shard"] != second_alt["ref_shard"]:  # Move in same stream
        # Move to upstream
        before_first_seq = before_first
        before_second_seq = before_second[:-1]  # Exclude ref
        if first_alt["ref_shard"] == "down":
            before_first_seq = before_first[:-1]  # Exclude ref
            before_second_seq = before_second
        for nt_first, nt_second in zip(before_first_seq[::-1], before_second_seq[::-1]):
            if nt_first != nt_second:
                break
            cipos_start -= 1
        # Move to downstream
        after_first_seq = after_first[1:]  # Exclude ref
        after_second_seq = after_second
        if first_alt["ref_shard"] == "down":
            after_first_seq = after_first
            after_second_seq = after_second[1:]  # Exclude ref
        for nt_first, nt_second in zip(after_first_seq, after_second_seq):
            if nt_first != nt_second:
                break
            cipos_end += 1
        # Update records
        if cipos_start != 0 or cipos_end != 0:
            first.pos = first.pos + cipos_start
            first.ref = seq_handler.getSub(first.chrom, first.pos, first.pos)
            first.info["CIPOS"] = [0, cipos_end - cipos_start]
            second.pos = first_alt["pos"] + cipos_start
            second.ref = seq_handler.getSub(second.chrom, second.pos, second.pos)
            second.info["CIPOS"] = first.info["CIPOS"]
            first.alt[0] = encodedAlt(
                first_alt | {"pos": second.pos, "nt": first.ref}
            )
            second.alt[0] = encodedAlt(
                second_alt | {"pos": first.pos, "nt": second.ref}
            )
    else:  # Move in opposite stream (same shard)
        # Move before first cointaining breakend and after second excluding breakend
        before_first_seq = before_first  # Contain ref
        after_second_seq = after_second[1:]  # Exclude ref
        if first_alt["ref_shard"] == "down":
            before_first_seq = before_first[:-1]  # Exclude ref
            after_second_seq = after_second  # Contain ref
        for nt_first, nt_second in zip(before_first_seq[::-1], after_second_seq):
            if nt_first != nt_second:
                break
            cipos_start -= 1
        # Move before second cointaining breakend and after first excluding breakend
        after_first_seq = after_first[1:]  # Exclude ref
        before_second_seq = before_second  # Contain ref
        if first_alt["ref_shard"] == "down":
            after_first_seq = after_first  # Contain ref
            before_second_seq = before_second[:-1]  # Exclude ref
        for nt_first, nt_second in zip(after_first_seq, before_second_seq[::-1]):
            if nt_first != nt_second:
                break
            cipos_end += 1
        # Update records
        if cipos_start != 0 or cipos_end != 0:
            first.pos = first.pos + cipos_start
            first.ref = seq_handler.getSub(first.chrom, first.pos, first.pos)
            first.info["CIPOS"] = [0, cipos_end - cipos_start]
            second.pos = first_alt["pos"] - cipos_end  # because cipos_start for first is - cipos_end for second
            second.ref = seq_handler.getSub(second.chrom, second.pos, second.pos)
            second.info["CIPOS"] = first.info["CIPOS"]
            second_down_pos = second.pos + (cipos_end - cipos_start)
            first_down_pos = first.pos + (cipos_end - cipos_start)
            first.alt[0] = encodedAlt(
                first_alt | {"pos": second_down_pos, "nt": first.ref}
            )
            second.alt[0] = encodedAlt(
                second_alt | {"pos": first_down_pos, "nt": second.ref}
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Replace N in alt and ref by the convenient nucleotid and move each breakend in pair at the left most position and add uncertainty iin CIPOS tag.')
    parser.add_argument('-p', '--sequence-padding', type=int, default=100, help='Number of nucleotids to inspect before and after the breakends: upstream and downstream movements are limited to this number of nucleotids. [Default: %(default)s]')
    parser.add_argument('-t', '--trace-unstandard', action='store_true', help='Use this option to add "UNSTD" tag in record INFO. This tag contains the trace of the variant before standardization: chromosome:position=reference/alternative.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the fusions file (format: VCF).')
    group_input.add_argument('-g', '--input-genome', required=True, help='Genome reference used in fusion calling to produced the inputed VCF (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the outputted file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    with IdxFastaIO(args.input_genome) as genome_reader:
        with BreakendVCFIO(args.input_variants) as reader:
            with BreakendVCFIO(args.output_variants, "w") as writer:
                writer.copyHeader(reader)
                writer.info["CIPOS"] = HeaderInfoAttr("CIPOS", type="Integer", number="2", description="Confidence interval around POS")
                if args.trace_unstandard:
                    writer.info["UNSTD"] = HeaderInfoAttr("UNSTD", type="String", number="1", description="Breakend id (chromosome:position=reference/alternative) before standardization")
                writer.writeHeader()
                for first, second in reader:
                    if len(first.alt) != 1 or len(second.alt) != 1:
                        raise Exception(f"Standardiztion cannot be used on multi-allelic variant {first.getName()}.")
                    if "IMPRECISE" not in set(first.info) and len(decodedAlt(first.alt[0])["nt"]) == 1:
                        if args.trace_unstandard:
                            first.info["UNSTD"] = first.getName()
                            second.info["UNSTD"] = second.getName()
                        fastStandardize(first, second, genome_reader, args.sequence_padding)
                    replaceUnknownNt(first, genome_reader)
                    replaceUnknownNt(second, genome_reader)
                    writer.write(first, second)
    log.info("End of job")
