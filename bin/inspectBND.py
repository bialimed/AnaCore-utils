#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.bed import getSortedAreasByChr
from anacore.fusion import BreakendVCFIO
from anacore.gff import GFF3IO
from anacore.gtf import loadModel
from anacore.region import iterOverlapped, RegionList
import argparse
import json
import logging
import numpy as np
import os
import pysam
from statistics import mean, median
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def annotateDomains(annot_by_gene_id, in_domains):
    """
    Add proteins domains information in annotations.

    :param annot_by_gene_id: Light genomic annotations by gene ID. Only genes overlapping breakends are present.
    :type annot_by_gene_id: dict
    :param in_domains: Path to the domains annotations file (format: GFF3). Each entry in file is one proteic domain. Required attributes are "Note" (long name) and "target_protein" (ID of the protein present in "input-annotations"). See AnaCore-utils/bin/ensemblInterProToGFF.py.
    :type in_domains: str
    """
    # Get protein by ID
    prot_by_id = dict()
    for gene_id, gene in annot_by_gene_id.items():
        for transcript in gene["transcripts"]:
            for protein in transcript["proteins"]:
                prot_by_id[protein["name"]] = protein
    # Add domains
    with GFF3IO(in_domains) as reader_dom:
        for domain in reader_dom:
            prot_id = domain.annot["target_protein"].split(".", 1)[0]
            if prot_id in prot_by_id:
                protein = prot_by_id[prot_id]
                if "domains" not in protein["annot"]:
                    protein["annot"]["domains"] = []
                domain_data = {
                    "annot": {"desc": domain.annot["Note"]},
                    "end": domain.end,
                    "start": domain.start
                }
                if "Dbxref" in domain.annot:
                    domain_data["annot"]["id"] = domain.annot["Dbxref"]
                if "sub_segments" in domain.annot:
                    domain_data["annot"]["sub"] = domain.annot["sub_segments"]
                protein["annot"]["domains"].append(domain_data)


def filteredByOverlap(targets_by_chr, selected_genes):
    """
    Return targeted areas overlapping selected genes.

    :param targets_by_chr: RegionList of targets by chromosomes.
    :type targets_by_chr: dict
    :param selected_genes: Translocated genes. Each element is a Gene object from genomic annotations.
    :type selected_genes: list
    :return: Targeted areas overlapping selected genes.
    :rtype: dict
    """
    # Genes by chromosome
    genes_by_chr = dict()
    for gene in selected_genes:
        chrom = gene.reference.name
        if chrom in selected_genes:
            genes_by_chr[chrom].append(gene)
        else:
            genes_by_chr[chrom] = RegionList([gene])
    # Find overlaps between targets and selected genes
    trimmed_targets_by_chr = dict()
    for chrom, genes in genes_by_chr.items():
        overlaps = list()
        for gene, targets in iterOverlapped(genes_by_chr[chrom], targets_by_chr[chrom]):
            for curr in targets:
                overlaps.append([
                    max(gene.start, curr.start),
                    min(gene.end, curr.end)
                ])
        consolidated_overlaps = list()
        if len(overlaps) > 0:
            overlaps = sorted(overlaps, key=lambda x: (x[0], x[1]))
            consolidated_overlaps = [overlaps[0]]
            prev = overlaps[0]
            for curr in overlaps[1:]:
                if curr[0] > prev[1]:
                    consolidated_overlaps.append(curr)
                    prev = curr
                else:
                    prev[1] = max(curr[1], prev[1])
        trimmed_targets_by_chr[chrom] = consolidated_overlaps
    return trimmed_targets_by_chr


def getAnnotByGene(tr_by_id, in_alignments, in_variants, min_base_qual=13, stranded=None, annotation_field="ANN"):
    """
    Return a light genomic annotations by gene ID of genes overlapping breakends.

    :param tr_by_id: Transcript object by ID from genomic annotations.
    :type tr_by_id: dict
    :param in_alignments: Path to the alignments file (format: BAM).
    :type in_alignments: str
    :param in_variants: Path to the annotated fusions file (format: VCF).
    :type in_variants: str
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param stranded: If your RNA library preps are stranded choose the read matching the transcript strand (R1: sense prep, R2: antisense prep).
    :type stranded: None|boolean
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    :return: A light genomic annotations by gene ID of genes overlapping breakends.
    :rtype: dict
    """
    processed_tr = set()
    annot_by_gene_id = {}
    with BreakendVCFIO(in_variants, annot_field=annotation_field) as reader:
        for first, second in reader:
            for bnd in [first, second]:
                for curr_annot in bnd.info[annotation_field]:
                    gene_id = curr_annot["Gene"]
                    tr_id = curr_annot["Feature"]
                    if tr_id not in processed_tr:
                        processed_tr.add(tr_id)
                        tr = tr_by_id[tr_id]
                        exons = []
                        for curr_exon in tr.children:
                            annotations = {}
                            if stranded is None:
                                depths = getDepths(in_alignments, curr_exon, min_base_qual)
                                annotations["depth"] = {
                                    "med": median(depths),
                                    "mean": mean(depths)
                                }
                            else:
                                depths_R1_fwd, depths_R2_fwd = getStrandedDepths(in_alignments, curr_exon, min_base_qual)
                                annotations["depth"] = {
                                    "med": [median(depths_R1_fwd), median(depths_R2_fwd)],
                                    "mean": [mean(depths_R1_fwd), mean(depths_R2_fwd)]
                                }
                                if stranded == "R2":
                                    annotations["depth"]["med"].reverse()
                                    annotations["depth"]["mean"].reverse()
                            exons.append({
                                "start": curr_exon.start,
                                "end": curr_exon.end,
                                "annot": annotations
                            })
                        # Add to gene
                        if gene_id not in annot_by_gene_id:
                            annot_by_gene_id[gene_id] = {
                                "name": gene_id,
                                "transcripts": []
                            }
                        annot_by_gene_id[gene_id]["transcripts"].append({
                            "name": tr_id,
                            "strand": tr.strand,
                            "proteins": [{"name": curr_prot.annot["id"], "start": curr_prot.start, "end": curr_prot.end, "annot": {}} for curr_prot in tr.proteins],
                            "exons": exons
                        })
    return annot_by_gene_id


def getDepths(in_aln, region, min_base_qual=13, excluded_flags="UNMAP,SECONDARY,QCFAIL,DUP"):
    """
    Return list of depths on region.

    :param in_aln: Path to the alignments file.
    :type in_aln: str
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :param excluded_flags: Discard any read that has any of the flags specified in the comma-separated list.
    :type excluded_flags: str
    :return: Depths on region.
    :rtype: list
    """
    out_str = pysam.depth(
        "-a",
        "-G", excluded_flags,
        "-J",
        "-q", str(min_base_qual),
        "-r", "{}:{}-{}".format(region.reference.name, region.start, region.end),
        in_aln
    )
    return [int(elt.split()[2]) for elt in out_str.strip().split("\n")]


def getStrandedDepths(in_aln, region, min_base_qual=13):
    """
    Return lists of stranded depths on region. The first is list of depths from R1 forward. The second is list of depths from R1 reverse.

    :param in_aln: Path to the alignments file.
    :type in_aln: str
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :param min_base_qual: Minimum base quality to count this read base on depth.
    :type min_base_qual: int
    :return: Depths on region.
    :rtype: list
    """
    excluded_flags = ["UNMAP", "SECONDARY", "QCFAIL", "DUP"]
    # R1
    r1_total = getDepths(
        in_aln, region, min_base_qual,
        ",".join(excluded_flags + ["READ2"])
    )
    r1_forward = getDepths(
        in_aln, region, min_base_qual,
        ",".join(excluded_flags + ["READ2", "REVERSE"])
    )
    r1_reverse = np.subtract(r1_total, r1_forward)
    # R2
    r2_total = getDepths(
        in_aln, region, min_base_qual,
        ",".join(excluded_flags + ["READ1"])
    )
    r2_forward = getDepths(
        in_aln, region, min_base_qual,
        ",".join(excluded_flags + ["READ1", "REVERSE"])
    )
    r2_reverse = np.subtract(r2_total, r2_forward)
    return (
        list(map(int, np.add(r1_forward, r2_reverse))),
        list(map(int, np.add(r1_reverse, r2_forward)))
    )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Produce data to inspect fusions breakends.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
    parser.add_argument('-q', '--min-base-qual', default=10, help="Minimum quality to take a read base into account in depth calculation. [Default: %(default)s]")
    parser.add_argument('-s', '--stranded', choices=["R1", "R2"], help='If your RNA library preps are stranded choose the read matching the transcript strand (R1: sense prep, R2: antisense prep).')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-l', '--input-alignments', required=True, help='Path to the alignments file (format: BAM).')
    group_input.add_argument('-a', '--input-annotations', required=True, help='Path to the genomic annotations file (format: GTF). Genes, transcripts and proteins used in fusion annotation process.')
    group_input.add_argument('-d', '--input-domains', help='Path to the domains annotations file (format: GFF3). Each entry in file is one proteic domain. Required attributes are "Note" (long name) and "target_protein" (ID of the protein present in "input-annotations"). See AnaCore-utils/bin/ensemblInterProToGFF.py.')
    group_input.add_argument('-t', '--input-targets', help='Path to the targets file (format: BED). If your protocol is targeted.')
    group_input.add_argument('-i', '--input-variants', required=True, help='Path to the annotated fusions file (format: VCF).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-annotations', required=True, help='Path to fusions details (format: JSON).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load annotations
    log.info("Load model from {}.".format(args.input_annotations))
    transcripts = loadModel(args.input_annotations, "transcripts")
    tr_by_id = {curr_tr.annot["id"].split(".")[0]: curr_tr for curr_tr in transcripts}

    # Get annotations by gene
    log.info("Get fusions genomic annotations from {}.".format(args.input_variants))
    annot_by_gene_id = getAnnotByGene(
        tr_by_id, args.input_alignments, args.input_variants, args.min_base_qual,
        args.stranded, args.annotation_field
    )

    # Annotate domains area on proteins
    if args.input_domains is not None:
        log.info("Annotate proteins domains areas on proteins from {}.".format(args.input_domains))
        annotateDomains(annot_by_gene_id, args.input_domains)

    # Store targeted areas overlapping selected genes
    selected_targets_by_chr = None
    if args.input_targets is not None:
        log.info("Store selected targets from {}.".format(args.input_targets))
        targets_by_chr = getSortedAreasByChr(args.input_targets)
        selected_genes = [curr_tr.parent for curr_tr in transcripts if curr_tr.parent.annot["id"].split(".")[0] in annot_by_gene_id]
        selected_targets_by_chr = filteredByOverlap(targets_by_chr, selected_genes)

    # Write output
    log.info("Write output.")
    with open(args.output_annotations, "w") as writer:
        json.dump({"gene_by_id": annot_by_gene_id, "targets_by_chr": selected_targets_by_chr}, writer)
    log.info("End of job")
