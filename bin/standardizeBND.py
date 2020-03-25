#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import re
import sys
import logging
import argparse
from anacore.vcf import HeaderInfoAttr
from anacore.fusion import BreakendVCFIO, getAltFromCoord, getCoordStr, getStrand
from anacore.sequenceIO import IdxFastaIO


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
    if "N" in record.alt[0]:
        matches = re.search("[\[\]](.+):(\d+)[\[\]]", record.alt[0])
        mate_chrom = matches.group(1)
        mate_pos = int(matches.group(2))
        mate_nt = seq_handler.getSub(mate_chrom, mate_pos, mate_pos)
        record.alt[0] = record.alt[0].replace("N", mate_nt)


def fastStandardize(first, second, seq_handler, padding=100):
    """
    Each breakend of the pair is placed at the left most position, and the uncertainty is represented with the CIPOS tag. The ALT string is then constructed assuming this choice.

    :param first: The breakend of the first shard in fusion (donor).
    :type first: anacore.vcf.VCFRecord
    :param second: The breakend of the second shard in fusion (acceptor).
    :type second: anacore.vcf.VCFRecord
    :param seq_handler: Indexed reader for the reference genome used in fusion calling.
    :type seq_handler: anacore.sequenceIO.IdxFastaIO
    :param padding: *******************************
    :type padding: int
    """
    first_strand = getStrand(first, True)
    second_strand = getStrand(second, False)
    before_first = seq_handler.getSub(first.chrom, max(first.pos - padding, 1), first.pos)
    before_second = seq_handler.getSub(second.chrom, max(second.pos - padding, 1), second.pos)
    after_first = seq_handler.getSub(first.chrom, first.pos, first.pos + padding)
    after_second = seq_handler.getSub(second.chrom, second.pos, second.pos + padding)
    cipos_start = 0
    cipos_end = 0
    if first_strand == second_strand:  # Same strand
        # Move to upstream
        before_first_seq = before_first
        before_second_seq = before_second[:-1]
        if first_strand == "-":
            before_first_seq = before_first[:-1]
            before_second_seq = before_second
        for nt_first, nt_second in zip(before_first_seq[::-1], before_second_seq[::-1]):
            if nt_first != nt_second:
                break
            cipos_start -= 1
        # Move to downstream
        after_first_seq = after_first[1:]
        after_second_seq = after_second
        if first_strand == "-":
            after_first_seq = after_first
            after_second_seq = after_second[1:]
        for nt_first, nt_second in zip(after_first_seq, after_second_seq):
            if nt_first != nt_second:
                break
            cipos_end += 1
        # Update records
        if cipos_start != 0 or cipos_end != 0:
            first.pos = first.pos + cipos_start
            first.ref = seq_handler.getSub(first.chrom, first.pos, first.pos)
            first.info["CIPOS"] = [0, cipos_end - cipos_start]
            second.pos = second.pos + cipos_start
            second.ref = seq_handler.getSub(second.chrom, second.pos, second.pos)
            second.info["CIPOS"] = [0, cipos_end - cipos_start]
            first.alt[0], second.alt[0] = getAltFromCoord(
                getCoordStr(first, True),
                getCoordStr(second, False)
            )
    else:  # Different strand
        before_second = getComplement(before_second)
        after_second = getComplement(after_second)
        # Move before first cointaining breakend and after second excluding breakend
        before_first_seq = before_first
        after_second_seq = after_second[1:]
        for nt_first, nt_second in zip(before_first_seq[::-1], after_second_seq):
            if nt_first != nt_second:
                break
            cipos_start -= 1
        # Move before second cointaining breakend and after first excluding breakend
        after_first_seq = after_first[1:]
        before_second_seq = before_second
        for nt_first, nt_second in zip(after_first_seq, before_second_seq[::-1]):
            if nt_first != nt_second:
                break
            cipos_end += 1
        # Update records
        if cipos_start != 0 or cipos_end != 0:
            first.pos = first.pos + cipos_start
            first.ref = seq_handler.getSub(first.chrom, first.pos, first.pos)
            first.info["CIPOS"] = [0, cipos_end - cipos_start]
            second.pos = second.pos - cipos_end  # because cipos_start for first is - cipos_end for second
            second.ref = seq_handler.getSub(second.chrom, second.pos, second.pos)
            second.info["CIPOS"] = first.info["CIPOS"]
            second_down_pos = second.pos + (cipos_end - cipos_start)
            first.alt[0], trash = getAltFromCoord(
                getCoordStr(first, True),
                {"chrom": second.chrom, "pos": second_down_pos, "strand": second_strand}
            )
            first_down_pos = first.pos + (cipos_end - cipos_start)
            trach, second.alt[0] = getAltFromCoord(
                {"chrom": first.chrom, "pos": first_down_pos, "strand": first_strand},
                getCoordStr(second, False)
            )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Replace N in alt and ref by the convenient nucleotid and move each breakend in pair at the left most position and add uncertainty iin CIPOS tag.')
    parser.add_argument('-t', '--trace-unstandard', action='store_true', help='Use this option to add "UNSTD" tag in record INFO. This tag contains the trace of the variant before standardization: chromosome:position=reference/alternative.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the fusions file (format: VCF).')
    group_input.add_argument('-g', '--input-genome', required=True, help='Genome reference used in fusion calling to produced the inputed VCF (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
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
                    if args.trace_unstandard:
                        first.info["UNSTD"] = "{}:{}={}/{}".format(first.chrom, first.pos, first.ref, "/".join(first.alt))
                        second.info["UNSTD"] = "{}:{}={}/{}".format(second.chrom, second.pos, second.ref, "/".join(second.alt))
                    if "Imprecise" not in set(first.filter):
                        fastStandardize(first, second, genome_reader)
                    replaceUnknownNt(first, genome_reader)
                    replaceUnknownNt(second, genome_reader)
                    writer.write(first, second)
    log.info("End of job")
