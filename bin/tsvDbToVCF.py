#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2024 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'


from anacore.db.homo_sapiens.accession import ChrAccession
from anacore.sequenceIO import IdxFastaIO
from anacore.sv import HashedSVIO
from anacore.vcf import HeaderInfoAttr, VCFIO, VCFRecord
import argparse
import logging
import os
import pysam
import pysam.bcftools
import sys
import yaml


################################################################################
#
# PRESETS
#
################################################################################
ALPHAMISSENSE = '''variant:
    region: "#CHROM"
    start: POS
    ref: REF
    alt: ALT
    id:
info:
    am_pathogenicity:
        name: score
        description: "Calibrated AlphaMissense pathogenicity scores (ranging between 0 and 1), which can
be interpreted as the predicted probability of a variant being clinically pathogenic."
        type: Float
    am_class:
        name: class
        description: "Pathogenicity class. Three discrete categories: likely_benign, ambiguous or likely_pathogenic. These are derived using the following thresholds of am_pathogenicity: likely_benign if am_pathogenicity < 0.34; likely_pathogenic if am_pathogenicity > 0.564; ambiguous otherwise."'''

CADD = '''variant:
    region: "#Chrom"
    start: Pos
    ref: Ref
    alt: Alt
    id:
info:
    RawScore:
        description: Pathogenicity prediction score.
        type: Float
    PHRED:
        description: Pathogenicity prediction score after normalisation and PHRED scaling.
        type: Float'''

DBNSFP = '''variant:
    region: "#chr"
    start: pos(1-based)
    ref: ref
    alt: alt
    id: rs_dbSNP
info:
    AlphaMissense_pred:
        name: AlphaMissense_class
        description: Pathogenicity prediction class.
    AlphaMissense_rankscore:
        description: Pathogenicity prediction score after scaling.
    CADD_phred:
        description: Pathogenicity prediction score after normalisation and PHRED scaling.
        type: Float
    CADD_raw_rankscore:
        name: CADD_rankscore
        description: Pathogenicity prediction score after scaling.
        type: Float
    clinvar_clnsig:
        name: ClinVar_class
        description: Pathogenicity prediction class.
        type: Float
    MetaLR_rankscore:
        description: Pathogenicity prediction score after scaling.
        type: Float
    MetaSVM_rankscore:
        description: Pathogenicity prediction score after scaling.
        type: Float
    VEST4_rankscore:
        description: Pathogenicity prediction score after scaling.
        type: Float'''


################################################################################
#
# FUNCTIONS
#
################################################################################
def convert(cfg, args):
    """
    Convert TSV variants database to VCF.

    :param cfg: Convert configuration. Database must be describe in two keys. First: "variant" contains column titles to find chrom, pos, ref, alt and optionaly id. Second: "info" contains the list of columns with override behaviour from default.
    :type cfg: dict
    :param log: The logger instance.
    :type log: logging.Logger
    """
    with IdxFastaIO(args.input_sequences) as seq_reader, HashedSVIO(args.input_variants) as reader, VCFIO(args.output_variants + ".tmp", "w") as writer:
        if args.convert_chrom:
            new_chrom_by_chrom = getRenameRules(seq_reader.index)
        # Header
        if args.source:
            writer.extra_header.append("##source={}".format(args.source))
        for chr_name, chr_info in sorted(seq_reader.index.items()):
            writer.extra_header.append(
                "##contig=<ID={},length={}>".format(chr_name, chr_info.length)
            )
        if args.discard_undef:
            selected_col = cfg.get("info", {}).keys()
        else:
            selected_col = set(reader.titles) - set(cfg["variant"].values())
        vcf_info_fields = dict()
        for tsv_name in selected_col:
            if tsv_name not in cfg.get("info", {}):  # Not defined field but keep undef
                writer.info[tsv_name] = HeaderInfoAttr(tsv_name, "", "String", "1")
                vcf_info_fields[tsv_name] = tsv_name
            else:  # Defined field
                field_cfg = cfg["info"][tsv_name]
                vcf_name = field_cfg.get("name", tsv_name)
                writer.info[vcf_name] = HeaderInfoAttr(
                    vcf_name, field_cfg.get("description", ""),
                    field_cfg.get("type", "String"), field_cfg.get("number", "1")
                )
                vcf_info_fields[tsv_name] = vcf_name
        writer.writeHeader()
        # Records
        for rec in reader:
            if "," in rec[cfg["variant"]["alt"]]:
                raise Exception(
                    "Multiple alternative variants are not allowed: {} in line {}".format(
                        rec[cfg["variant"]["alt"]], reader.current_line_nb
                    )
                )
            chrom = rec[cfg["variant"]["region"]]
            if args.convert_chrom:
                chrom = new_chrom_by_chrom[ChrAccession.toHumanName(chrom)]
            vcf_rec = VCFRecord(
                region=chrom,
                position=int(rec[cfg["variant"]["start"]]) + (1 if args.start_0_based else 0),
                refAllele=rec[cfg["variant"]["ref"]],
                altAlleles=[rec[cfg["variant"]["alt"]]]
            )
            if cfg["variant"].get("id"):
                vcf_rec.id = rec[cfg["variant"]["id"]]
            for tsv_name, vcf_name in vcf_info_fields.items():
                if rec[tsv_name].strip() != "":
                    vcf_rec.info[vcf_name] = rec[tsv_name]
            vcf_rec.standardize(seq_reader, 500)
            writer.write(vcf_rec)


def getCfgFromPreset(preset):
    """
    Return configuration from selected preset.

    :param preset: Selected preset.
    :type preset: str
    :return: Configuration from selected preset.
    :rtype: dict
    """
    preset_content = None
    if preset == "AlphaMissense":
        preset_content = ALPHAMISSENSE
    elif preset == "CADD":
        preset_content = CADD
    elif preset == "dbNSFP":
        preset_content = DBNSFP
    else:
        raise ValueError("Preset {} is not valid.".format(preset))
    return yaml.safe_load(preset_content)


def getRenameRules(seq_chroms):
    """
    Return for human genome the match between the chromosome name in database (value) and the chromosome name in sequences file (key).

    :Example:
        If sequences file contains:
            * **refseq** ID like NC_012920.9, the function return **{"MT": "NC_012920.9", ...}**
            * **genbank** ID like J01415, the function return **{"MT": "NC_012920.9", ...}**
            * **chr prefixed** ID like chrMT, the function return **{"MT": "chrMT", ...}**
            * **short mitochondiral chromosome notation** M, the function return **{"MT": "M", ...}**

    :param seq_chroms: List of chromosomes name in sequences file.
    :type seq_chroms: iterable
    :return: Match between the chromosome name in database (value) and the chromosome name in sequences file (key).
    :rtype: dict
    """
    rename_by_chrom = dict()
    for new_name in seq_chroms:
        old_name = new_name
        if old_name.startswith("chr"):
            old_name = new_name[3:]
        if old_name == "M":
            old_name = "MT"
        rename_by_chrom[old_name] = new_name
    return rename_by_chrom


################################################################################
#
# MAIN
#
################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Convert TSV variants database to indexed VCF.')
    parser.add_argument('-0', '--start-0-based', action='store_true', help='Start position in input-variants is 0-based.')
    parser.add_argument('-c', '--convert-chrom', action='store_true', help='Convert human chromosome accession to human readable name used in sequences file.')
    parser.add_argument('-d', '--discard-undef', action='store_true', help='Do not keep columns not defined in input-configuration/preset. Otherwise, all columns from TSV are add as INFO in VCF.')
    parser.add_argument('-u', '--source', help='Add source tag in meta-information of VCF.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_config = parser.add_mutually_exclusive_group(required=True)
    group_config.add_argument('-f', '--input-configuration', help='Path to configuration specify fields for variant coordinates and INFO tags (format: YAML). Database must be describe in two parts. First: "variant" contains column titles to find chrom, pos, ref, alt and optionaly id. Second: "info" contains the list of columns with override behaviour from default. Default behaviour for each column is to convert in INFO tag with name corresponding to column name, blank description, type=String and number=1. In info section each entry use TSV column tilte as key and can set name (INFO tag name), description, type (Integer, Float, Character and String) and number (1, A, G, R, and .).')
    group_config.add_argument('-p', '--preset', choices=["AlphaMissense", "CADD", "dbNSFP"], help='Preset for converting configuration.')
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to variants file (format: TSV).')
    group_input.add_argument('-s', '--input-sequences', required=True, help='Path to reference sequences file (format: fasta with faidx).')
    group_output = parser.add_argument_group('Outputs')
    group_input.add_argument('-o', '--output-variants', required=True, help='Path to variants file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Convert
    cfg = dict()
    if args.input_configuration:
        with open(args.input_configuration) as reader:
            cfg = yaml.safe_load(reader)
    else:
        cfg = getCfgFromPreset(args.preset)
    convert(cfg, args)

    # Index
    pysam.bcftools.sort('-o', args.output_variants + ".tmp_sort", args.output_variants + ".tmp", catch_stdout=False)
    os.remove(args.output_variants + ".tmp")
    pysam.tabix_compress(args.output_variants + ".tmp_sort", args.output_variants)
    os.remove(args.output_variants + ".tmp_sort")
    pysam.tabix_index(args.output_variants, preset="vcf")
    log.info("End of job")
