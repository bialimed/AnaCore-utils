#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
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


def getAnnotByGene(tr_by_id, in_alignments, in_variants, stranded=None, annotation_field="ANN"):
    """
    Return a light genomic annotations by gene ID of genes overlapping breakends.

    :param tr_by_id: Transcript object by ID from genomic annotations.
    :type tr_by_id: dict
    :param in_alignments: Path to the alignments file (format: BAM).
    :type in_alignments: str
    :param in_variants: Path to the annotated fusions file (format: VCF).
    :type in_variants: str
    :param stranded: If your RNA library preps are stranded choose the read matching the transcript strand (R1: sense prep, R2: antisense prep).
    :type stranded: None|boolean
    :param annotation_field: Field used for store annotations.
    :type annotation_field: str
    :return: A light genomic annotations by gene ID of genes overlapping breakends.
    :rtype: dict
    """
    processed_tr = set()
    annot_by_gene_id = {}
    with pysam.AlignmentFile(in_alignments, "rb") as reader_aln:
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
                                    depths = getDepths(reader_aln, curr_exon)
                                    annotations["depth"] = {
                                        "med": median(depths),
                                        "mean": mean(depths)
                                    }
                                else:
                                    depths_R1_fwd, depths_R2_fwd = getStrandedDepths(reader_aln, curr_exon)
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
                                    "transcripts": [],
                                    "annot": {}
                                }
                            annot_by_gene_id[gene_id]["transcripts"].append({
                                "name": tr_id,
                                "strand": tr.strand,
                                "proteins": [{"name": curr_prot.annot["id"], "start": curr_prot.start, "end": curr_prot.end, "annot": {}} for curr_prot in tr.proteins],
                                "exons": exons
                            })
    return annot_by_gene_id


def getDepths(reader_aln, region):
    """
    Return list of depths on region.

    :param reader_aln: The file handle to the alignments file.
    :type reader_aln: pysam.AlignmentFile
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :return: Depths on region.
    :rtype: list
    """
    depths = []
    expected_pos = region.start
    for aln_col in reader_aln.pileup(
        region.reference.name,
        region.start - 1,  # 0-based
        region.end,  # 1-based
        truncate=True,
        max_depth=1000000
    ):  # By defaul filter out unmap, secondary, qcfail and duplicate
        pileup_pos = aln_col.pos + 1  # 0-based
        while expected_pos != pileup_pos:  # Take into account positions without alingment
            depths.append(0)
            expected_pos += 1
        curr_dp = sum([1 for read in aln_col.pileups])
        depths.append(curr_dp)
        expected_pos += 1
    end_limit = region.end + 1  # expected_pos is the current unprocessed position
    while expected_pos != end_limit:  # Take into account positions without alingment
        depths.append(0)
        expected_pos += 1
    return depths


def getStrandedDepths(reader_aln, region):
    """
    Return lists of stranded depths on region. The first is list of depths from R1 forward. The second is list of depths from R1 reverse.

    :param reader_aln: The file handle to the alignments file.
    :type reader_aln: pysam.AlignmentFile
    :param region: The evaluated region.
    :type region: anacore.region.Region
    :return: Stranded depths on region. The first is list of depths from R1 forward. The second is list of depths from R1 reverse.
    :rtype: (list, list)
    """
    depths_R1_foward = []
    depths_R1_reverse = []
    expected_pos = region.start
    for aln_col in reader_aln.pileup(
        region.reference.name,
        region.start - 1,  # 0-based
        region.end,  # 1-based
        truncate=True,
        max_depth=1000000
    ):  # By defaul filter out unmap, secondary, qcfail and duplicate
        pileup_pos = aln_col.pos + 1  # 0-based
        while expected_pos != pileup_pos:  # Take into account positions without alingment
            depths_R1_foward.append(0)
            depths_R1_reverse.append(0)
            expected_pos += 1
        curr_dp_R1_forward = 0
        curr_dp_R1_reverse = 0
        for read in aln_col.pileups:
            if read.is_read1:
                if read.is_reverse:
                    curr_dp_R1_reverse += 1
                else:
                    curr_dp_R1_forward += 1
            else:  # Read 2
                if read.is_reverse:
                    curr_dp_R1_forward += 1
                else:
                    curr_dp_R1_reverse += 1
        depths_R1_foward.append(curr_dp_R1_forward)
        depths_R1_reverse.append(curr_dp_R1_reverse)
        expected_pos += 1
    end_limit = region.end + 1  # expected_pos is the current unprocessed position
    while expected_pos != end_limit:  # Take into account positions without alingment
        depths_R1_foward.append(0)
        depths_R1_reverse.append(0)
        expected_pos += 1
    return depths_R1_foward, depths_R1_reverse


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Produce data to inspect fusions breakends.')
    parser.add_argument('-f', '--annotation-field', default="ANN", help='Field used for store annotations. [Default: %(default)s]')
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
    annot_by_gene_id = getAnnotByGene(tr_by_id, args.input_alignments, args.input_variants, args.stranded, args.annotation_field)

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
