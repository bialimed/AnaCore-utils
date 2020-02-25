#!/usr/bin/env python3

__author__ = 'Veronique Ivashchenko and Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import sys
import argparse
import logging
from anacore.sv import SVIO
from anacore.vcf import VCFIO


########################################################################
#
# FUNCTIONS
#
########################################################################
def uGetHeaderAttr(header_line):
    """
    Return an instance of HeaderAttr or a child corresponding to the VCF line.

    :param header_line: Declaration line declaration of an VCF attributes from the VCF header.
    :type header_line: str
    :return: An instance of HeaderAttr or a child corresponding to the VCF line.
    :rtype: HeaderAttr or a child
    """
    # Get content
    header_category, header_content = header_line.split("=", 1)  # ##INFO=<ID=AD,Version="1">
    header_category = header_category[2:]  # ##INFO to INFO
    header_content = header_content[1:-1]  # <ID=AD,Version="1"> to ID=AD,Version="1"
    # Get attributes
    attributes = {}
    stack = ""
    opened_quote = False
    for curr_char in header_content:
        if curr_char == "," and not opened_quote:
            key, val = stack.split("=", 1)
            attributes[key.lower()] = val
            stack = ""
            opened_quote = False
        elif curr_char == '"':
            if stack[-1] == '\\':
                stack = stack[:-1] + '"'  # replace '\"' by '"'
            else:
                if opened_quote:  # The quote is the second
                    opened_quote = False
                else:  # The quote is the first
                    opened_quote = True
        else:
            stack += curr_char
    if stack != "":
        key, val = stack.split("=", 1)
        attributes[key.lower()] = val
    # Return
    header_class_name = "HeaderAttr".format(header_category.capitalize())
    header_class = getattr(sys.modules["anacore.vcf"], header_class_name)
    return header_class(**attributes)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser('Convert chromosome specification in STAR_Fusion VCF file into standard chromosome specification.')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-variants', required=True, help="Path to STAR_Fusion VCF file.")
    group_input.add_argument('-n', '--input-names', help="Path to file containing rename rules (format: TSV). The first column contains the old name and second contains the new. Without this parameter, the prefix chr is removed from main chromosomes names of human assembly.")
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help="Path to STAR_Fusion_converted VCF output file.")
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Get naming rules
    new_names = {'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5', 'chr6': '6', 'chr7': '7', 'chr8': '8', 'chr9': '9', 'chr10': '10', 'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15', 'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20', 'chr21': '21', 'chr22': '22', 'chrX': 'X', 'chrY': 'Y', 'chrM': 'MT'}
    if args.input_names:
        with SVIO(args.input_names, "r", separator="\t", has_title=False) as reader:
            for record in reader:
                new_names[record[0]] = record[1]

    # Process
    with VCFIO(args.output_variants, "w") as writer:
        with VCFIO(args.input_variants, "r") as reader:
            # Header
            writer.copyHeader(reader)
            for idx, curr_header in enumerate(writer.extra_header):
                if curr_header.startswith("##contig"):
                    content = uGetHeaderAttr(curr_header)
                    old_id = content.id
                    if content.id in new_names:
                        new_id = new_names[old_id]
                        writer.extra_header[idx] = curr_header.replace(
                            "ID={},".format(old_id),
                            "ID={},".format(new_id)
                        )
            writer.writeHeader()
            # Variants
            for record in reader:
                if record.chrom in new_names:
                    record.chrom = new_names[record.chrom]
                if "SVTYPE" in record.info and record.info["SVTYPE"] == "BND":
                    for idx, alt in enumerate(record.alt):
                        alt_left = alt.split(":")[0]
                        if "[" in alt_left:
                            curr_chrom = alt_left.split("[")[1]
                            if curr_chrom in new_names:
                                new_chrom = new_names[curr_chrom]
                                record.alt[idx] = alt.replace("[" + curr_chrom + ":", "[" + new_chrom + ":")
                        elif "]" in alt_left:
                            curr_chrom = alt_left.split("]")[1]
                            if curr_chrom in new_names:
                                new_chrom = new_names[curr_chrom]
                                record.alt[idx] = alt.replace("]" + curr_chrom + ":", "]" + new_chrom + ":")
                writer.write(record)

    # Log process
    log.info("End of job")
