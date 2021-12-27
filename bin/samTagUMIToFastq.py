#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.sequenceIO import FastqIO, Sequence
import argparse
import logging
import os
from pysam import AlignmentFile
import sys


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Convert BAM with UMI in specific tag to fastq with UMI in reads ID (see Illumina's reads ID format).")
    parser.add_argument('-b', '--barcode-tag', default="BC", help='Tag used in alignment file to store the barcode. [Default: %(default)s]')
    parser.add_argument('-k', '--keep-qc-failed', action='store_true', help='Keep QC failed reads.')
    parser.add_argument('-q', '--qual-offset', default=33, type=int, help='Quality offset in reads. [Default: %(default)s]')
    parser.add_argument('-r', '--reads-barcode', help='Reads barcode.')
    parser.add_argument('-t', '--umi-qual-tag', default="QX", help='Tag used in alignment file to store the UMI quality. [Default: %(default)s]')
    parser.add_argument('-u', '--umi-tag', default="RX", help='Tag used in alignment file to store the UMI. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-aln', required=True, help='The path to the alignments file (format: BAM).')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-reads', required=True, help='The path to the outputted reads file (format: FASTQ).')
    group_output.add_argument('-2', '--output-reads-2', help='The path to the outputted reads file R2 (format: FASTQ).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    if args.output_reads_2:  # Write all reads in a pair of files (R1 and R2)
        with FastqIO(args.output_reads_2, "w") as writer_r2:
            with FastqIO(args.output_reads, "w") as writer_r1:
                with AlignmentFile(args.input_aln, "rb", check_sq=False) as reader:
                    for curr_read in reader.fetch(until_eof=True):
                        if not curr_read.is_secondary and not curr_read.is_supplementary:
                            if args.keep_qc_failed or not curr_read.is_qcfail:
                                barcode = args.reads_barcode
                                if barcode is None and curr_read.has_tag(args.barcode_tag):
                                    barcode = curr_read.get_tag(args.barcode_tag).replace("-", "+")
                                description = "{}:{}:0:{} {}={}".format(
                                    "1" if curr_read.is_read1 else "2",
                                    "Y" if curr_read.is_qcfail else "N",
                                    "" if barcode is None else barcode,
                                    args.umi_qual_tag,
                                    curr_read.get_tag(args.umi_qual_tag)
                                )
                                read = Sequence(
                                    curr_read.query_name + ":" + curr_read.get_tag(args.umi_tag),
                                    curr_read.get_forward_sequence(),
                                    description,
                                    "".join([chr(elt + args.qual_offset) for elt in curr_read.get_forward_qualities()])
                                )
                                if curr_read.is_read1:
                                    writer_r1.write(read)
                                else:
                                    writer_r2.write(read)
    else:  # Write all reads in a unique file
        with FastqIO(args.output_reads, "w") as writer:
            with AlignmentFile(args.input_aln, "rb", check_sq=False) as reader:
                for curr_read in reader.fetch(until_eof=True):
                    if not curr_read.is_secondary and not curr_read.is_supplementary:
                        if args.keep_qc_failed or not curr_read.is_qcfail:
                            barcode = args.reads_barcode
                            if barcode is None and curr_read.has_tag(args.barcode_tag):
                                barcode = curr_read.get_tag(args.barcode_tag).replace("-", "+")
                            description = "{}:{}:0:{} {}={}".format(
                                "1" if curr_read.is_read1 else "2",
                                "Y" if curr_read.is_qcfail else "N",
                                "" if barcode is None else barcode,
                                args.umi_qual_tag,
                                curr_read.get_tag(args.umi_qual_tag)
                            )
                            read = Sequence(
                                curr_read.query_name + ":" + curr_read.get_tag(args.umi_tag),
                                curr_read.get_forward_sequence(),
                                description,
                                "".join([chr(elt + args.qual_offset) for elt in curr_read.get_forward_qualities()])
                            )
                            writer.write(read)
    log.info("End of job")
