#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'


import os
import sys
import logging
import argparse
from anacore.fusion import BreakendVCFIO
from anacore.sv import HashedSVIO
from anacore.vcf import HeaderInfoAttr


########################################################################
#
# FUNCTIONS
#
########################################################################
def sourcesBySymbols(in_known):
    """
    Return sources descriptions by fusion ID from database.

    :param in_known: Path to the file containing known fusions (format: TSV). This file must contains 3 columns : 5prim_gene, 3_prim_gene and sources. 5prim_gene and 3prim_gene are symbol with the same master name of the name in GTF used for the annotation of breakends. sources is a string containing db1name:entryId,entryId|db2name:entryId (example: cosmic_91:1743,1745|chimerdb_pub-V4:3427,3428).
    :type in_known: str
    :return: sources descriptions (db1name:entryId,entryId|db2name:entryId) by fusion ID (5primSymbol_@_3primSymbol).
    :rtype: dict
    """
    sources_by_symbols = {}
    with HashedSVIO(in_known) as reader:
        for record in reader:
            fusion_id = "{}_@_{}".format(record["5prim_gene"], record["3prim_gene"])
            sources_by_symbols[fusion_id] = record["sources"]
    return sources_by_symbols


def annotated(first, second, sources_by_symbols, annotation_field):
    """
    Return known fusions with databases names and entries ID known.

    :param first: Breakend of the 5' shard of the fusion.
    :type first: anacore.vcf.VCFRecord
    :param second: Breakend of the 3' shard of the fusion.
    :type second: anacore.vcf.VCFRecord
    :param sources_by_symbols: sources descriptions (db1name:entryId,entryId|db2name:entryId) by fusion ID (5primSymbol_@_3primSymbol).
    :type sources_by_symbols: dict
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    :return: Known fusions with databases names and entries ID known.
    :rtype: list
    """
    known = []
    for first_gene in {elt["SYMBOL"] for elt in first.info[annotation_field]}:
        for second_gene in {elt["SYMBOL"] for elt in second.info[annotation_field]}:
            if first_gene != "" and second_gene != "":
                fusion_id = "{}_@_{}".format(first_gene, second_gene)
                if fusion_id in sources_by_symbols:  # Only partners and order in fusion is evaluated (strand is skipped)
                    known.append(fusion_id + "=" + sources_by_symbols[fusion_id])
    return known


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Annotate fusions known in fusions database with databases names and entries ID known.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-p', '--input-known-partners', required=True, help='Path to the file containing known fusions (format: TSV). This file must contains 3 columns : 5prim_gene, 3_prim_gene and sources. 5prim_gene and 3prim_gene are symbol with the same master name of the name in GTF used for the annotation of breakends. sources is a string containing db1name:entryId,entryId|db2name:entryId (example: cosmic_91:1743,1745|chimerdb_pub-V4:3427,3428).')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the file containing variants annotated (format: VCF). The process use only the SYMBOL field of each breakend.')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-variants', required=True, help='Path to the annotated file. (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load knowns
    log.info("Load known partners from {}.".format(args.input_known_partners))
    sources_by_symbols = sourcesBySymbols(args.input_known_partners)

    # Annot variants
    log.info("Annotate known fusions partners.")
    with BreakendVCFIO(args.output_variants, "w", args.annotation_field) as writer:
        with BreakendVCFIO(args.input_variants, "r", args.annotation_field) as reader:
            # Header
            writer.copyHeader(reader)
            writer.info["known_partners"] = HeaderInfoAttr(
                id="known_partners",
                type="String",
                number=".",
                description="Database containing the fusion of these gene. Format: 5primSymbol_@_3primSymbol=db1name:entryId,entryId|db2name:entryId (example: BCR_@_ABL1=cosmic_91:1743,1745|chimerdb_pub-V4:3427,3428)"
            )
            writer.writeHeader()
            # Records
            for first, second in reader:
                if "RNA_FIRST" in first.info:  # Stranded
                    first.info["known_partners"] = annotated(first, second, sources_by_symbols, args.annotation_field)
                if "PRED_RNA_FIRST" in first.info or "PRED_RNA_FIRST" in second.info:  # Stranded by annot
                    if "PRED_RNA_FIRST" in second.info:
                        second, first = first, second
                    first.info["known_partners"] = annotated(first, second, sources_by_symbols, args.annotation_field)
                else:  # Unstranded
                    first.info["known_partners"] = annotated(first, second, sources_by_symbols, args.annotation_field) + annotated(second, first, sources_by_symbols, args.annotation_field)
                writer.write(first, second)
    log.info("End of job")
