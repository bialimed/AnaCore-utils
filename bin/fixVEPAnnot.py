#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import argparse
from anacore.annotVcf import AnnotVCFIO
from anacore.vcf import getAlleleRecord, VCFIO
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def changeCosmicAnnotations(record, annot_field, cosmic_reader):
    """
    Replace non-allele specific COSMIC annotations produced by VEP to allelle-specific annotations.

    .. warning::
    Alleles not annotated by VEP will remain unannotated despite data in databank.

    :param record: Annotated VCF record from VEP.
    :type record: anacore.vcf.VCFRecord
    :param annot_field: Field used to store annotations.
    :type annot_field: str
    :param cosmic_reader: File handler open on COSMIC databank with mode 'i'.
    :type cosmic_reader: anacore.vcf.VCFIO
    """
    annot_chr = record.chrom.upper()
    annot_name_prefix = "{}:{}={}/".format(
        annot_chr[3:] if annot_chr.startswith("CHR") else annot_chr,
        record.pos,
        record.ref.upper()
    )
    # Get overlapping COSMIC variants
    db_overlapping = [elt for elt in cosmic_reader.getSub(
        record.chrom[3:] if record.chrom.startswith("chr") else record.chrom,
        int(record.refStart()),
        int(record.refEnd() + 0.5)
    )]
    # Replace COSMIC annotations
    for annot in record.info[annot_field]:
        annot_name = annot_name_prefix + annot["Allele"].upper()
        new_existing = []
        # Remove old COSMIC annotations
        if annot["Existing_variation"] is not None:
            for curr_exist in annot["Existing_variation"].split("&"):
                if not curr_exist.startswith("COS"):
                    new_existing.append(curr_exist)
        # Add new COSMIC annotations
        db_ids = set()
        for db_record in db_overlapping:
            if len(db_record.alt) == 1:
                if annot_name == db_record.getName().upper():
                    db_ids = db_ids | set(db_record.id.split(";"))
            else:
                for alt_idx, db_alt in enumerate(db_record.alt):
                    db_alt_record = getAlleleRecord(cosmic_reader, db_record, alt_idx)
                    if annot_name == db_alt_record.getName().upper():
                        db_ids = db_ids | set(db_record.id.split(";"))
        new_existing += sorted(db_ids)
        # Change existing variants
        if len(new_existing) != 0:
            annot["Existing_variation"] = "&".join(new_existing)
        else:
            annot["Existing_variation"] = None


def getDatabankVersion(cosmic_reader):
    """
    Return COSMIC databank version from a file handler.

    :param cosmic_reader: File handler open on COSMIC databank.
    :type cosmic_reader: anacore.vcf.VCFIO
    :return: COSMIC databank version.
    :rtype: str
    """
    cosmic_version = None
    for curr_head in cosmic_reader.extra_header:
        if curr_head.startswith("##source="):
            cosmic_version = curr_head.lower().split("cosmicv")[1]
    return cosmic_version


def getVEPAlt(ref, alt):
    """
    Return the alternative allele in same format as annotation allele in VEP.

    :param ref: The reference allele.
    :type ref: str
    :param alt: The alternative allele.
    :type alt: str
    :return: The alternative allele in same format as annotation allele in VEP.
    :rtype: str
    """
    alleles = [ref] + alt
    # Replace empty marker by empty string
    for idx, cur_allele in enumerate(alleles):
        if cur_allele == "-":
            alleles[idx] = ""
    # Select shorter allele
    shorter_allele = alleles[0]
    for current_alt in alleles[1:]:
        if len(current_alt) < len(shorter_allele):
            shorter_allele = current_alt
    # Trim alleles
    trim = True
    while len(shorter_allele) != 0 and shorter_allele != "" and trim:
        for cur_allele in alleles:
            if len(cur_allele) == 0:
                trim = False
            elif cur_allele[0] != shorter_allele[0]:
                trim = False
        if trim:
            shorter_allele = shorter_allele[1:]
            for idx, cur_allele in enumerate(alleles):
                alleles[idx] = cur_allele[1:]
    # Replace empty by empty_marker
    for idx, cur_allele in enumerate(alleles):
        if cur_allele == "":
            alleles[idx] = "-"
    return alleles[1:]


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Reverse normalisation produced by VEP in allele annotation field.')
    parser.add_argument('-a', '--annotations-field', default="ANN", help='Field used to store annotations. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-c', '--input-cosmic', help='The path to the variants known in COSMIC (format: VCF with tbi). This option replace non-allele specific cosmic annotation produce by VEP to allelle-specific annotation. Ensembl is unfortunately not licensed to redistribute allele-specific data for cosmic.')
    group_input.add_argument('-i', '--input-variants', required=True, help='The path to the variants file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-variants', required=True, help='The path to the output file (format: VCF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    cosmic_reader = None
    with AnnotVCFIO(args.output_variants, "w", annot_field=args.annotations_field) as FH_out:
        with AnnotVCFIO(args.input_variants, annot_field=args.annotations_field) as FH_in:
            # Header
            FH_out.copyHeader(FH_in)
            if args.input_cosmic:
                cosmic_reader = VCFIO(args.input_cosmic, "i")
                cosmic_version = getDatabankVersion(cosmic_reader)
                FH_out.extra_header.append("##COSMIC={}".format(cosmic_version))
            FH_out.writeHeader()
            # Records
            for record in FH_in:
                # To upper
                record.ref = record.ref.upper()
                record.alt = [alt.upper() for alt in record.alt]
                for annot in record.info[FH_in.annot_field]:
                    annot["Allele"] = annot["Allele"].upper()
                # Change alternative representation
                for alt_idx, alt in enumerate(record.alt):
                    alt_record = getAlleleRecord(FH_in, record, alt_idx)
                    vep_alt = getVEPAlt(alt_record.ref, alt_record.alt)[0]
                    for idx_ann, annot in enumerate(alt_record.info[FH_in.annot_field]):
                        if annot["Allele"] == vep_alt:
                            annot["Allele"] = alt_record.alt[0]
                # Replace cosmic annotations
                if args.input_cosmic:
                    changeCosmicAnnotations(record, FH_in.annot_field, cosmic_reader)
                FH_out.write(record)
    log.info("End of job")
