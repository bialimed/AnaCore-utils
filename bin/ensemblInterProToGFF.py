#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

from anacore.gff import GFF3IO, GFF3Record
from anacore.gtf import loadModel
from anacore.sv import HashedSVIO
import argparse
import logging
import os
import sys


########################################################################
#
# FUNCTIONS
#
########################################################################
def hasOverlap(start, end, second_region):
    """
    Return True if the area defined by start and end has an overlap with the second region.

    :param start: Start of the first region.
    :type start: int
    :param end: End of the first region.
    :type end: int
    :param second_region: The second region.
    :type second_region: Region
    :return: True if the region has an overlap with the second region.
    :rtype: bool
    """
    return not start > second_region.end and not end < second_region.start


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Convert InterProScan annotations from ensembl's biomart in GFF.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-a', '--input-annotations', help='Path to the genomic annotations file (format: GTF). It contains coordinates for transcripts, exon and cds.')
    group_input.add_argument('-d', '--input-domains', help="Path to the domains annotations file (format: TSV). It contains InterPro domains extracted from ensembl's biomart (mandatory fields: 'Transcript stable ID version', 'Protein stable ID version', 'Interpro ID', 'Interpro Short Description', 'Interpro Description', 'Interpro start', 'Interpro end').")
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-annotations', required=True, help='Path to the domains annotations file (format: GFF).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Load annotations
    log.info("Load model from {}.".format(args.input_annotations))
    tr_by_id = {tr.annot["id"]: tr for tr in loadModel(args.input_annotations, "transcripts")}

    # Parse and convert domains data
    log.info("Parse and convert domains data from {}.".format(args.input_domains))
    domains_by_tr_id = dict()
    with HashedSVIO(args.input_domains) as reader:
        for record in reader:
            if record['Interpro ID'] != "":
                record['Interpro start'] = int(record['Interpro start'])
                record['Interpro end'] = int(record['Interpro end'])
                tr_id = record['Transcript stable ID version'].split(".", 1)[0]
                if tr_id not in tr_by_id:
                    log.warning("The transcript {} is missing in {}.".format(tr_id, args.input_annotations))
                else:
                    domain_id = record['Interpro ID']
                    # Get genomic coordinates
                    transcript = tr_by_id[tr_id]
                    protein = transcript.proteins[0]
                    if len(transcript.proteins) > 1:
                        msg = "The transcript {} is linked with several proteins {}.".format(tr_id, [prot.annot["id"] for prot in transcript.proteins])
                        log.error(msg)
                        raise Exception(msg)
                    try:
                        genomic_start = protein.getPosOnRef(record['Interpro start'], 1)
                    except Exception:
                        log.warning("Incomplete CDS for {} (transcript: {}).".format(protein, tr_id))
                        genomic_start = protein.start if transcript.strand == "+" else protein.end
                    try:
                        genomic_end = protein.getPosOnRef(record['Interpro end'], 3)
                    except Exception:
                        log.warning("Incomplete CDS for {} (transcript: {}).".format(protein, tr_id))
                        genomic_end = protein.end if transcript.strand == "+" else protein.start
                    if genomic_start > genomic_end:
                        aux = genomic_start
                        genomic_start = genomic_end
                        genomic_end = aux
                    # Add info to transcript
                    if tr_id not in domains_by_tr_id:
                        domains_by_tr_id[tr_id] = dict()
                    if domain_id not in domains_by_tr_id[tr_id] or not hasOverlap(genomic_start, genomic_end, domains_by_tr_id[tr_id][domain_id]):
                        domains_by_tr_id[tr_id][domain_id] = GFF3Record(
                            seq_id=transcript.reference.name,
                            source="InterProScan",
                            type="SO:0000839",  # cf. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md and ftp://ftp.geneontology.org/pub/go/doc/GO.xrf_abbs
                            start=genomic_start,
                            end=genomic_end,
                            strand=transcript.strand,
                            attributes={
                                "Dbxref": "InterPro:" + record['Interpro ID'],  # cf. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md and ftp://ftp.geneontology.org/pub/go/doc/GO.xrf_abbs
                                "Name": record['Interpro Short Description'],
                                "Note": record['Interpro Description'],
                                "target_protein": record['Protein stable ID version'],
                                "target_transcript": record['Transcript stable ID version']
                            }
                        )
                    else:  # Overlap
                        prev_annot = domains_by_tr_id[tr_id][domain_id]
                        log.debug(
                            "Merge domains {} {}-{} (aa: {}-{}) and {}-{} of transcript {}.".format(
                                record['Interpro ID'],
                                genomic_start,
                                genomic_end,
                                record['Interpro start'],
                                record['Interpro end'],
                                prev_annot.start,
                                prev_annot.end,
                                tr_id
                            )
                        )
                        prev_annot.start = min(genomic_start, prev_annot.start)
                        prev_annot.end = max(genomic_end, prev_annot.end)
    domains = [domain for tr_id, domain_by_dom_id in domains_by_tr_id.items() for domain_id, domain in domain_by_dom_id.items()]
    domains = sorted(domains, key=lambda x: (x.reference.name, x.start, x.end))

    # Add split info
    log.info("Add sub-segment without intron(s) for domains.")
    for curr_domain in domains:
        tr_id = curr_domain.annot["target_transcript"].split(".", 1)[0]
        transcript = tr_by_id[tr_id]
        segments = []
        for exon in sorted(transcript.children, key=lambda x: (x.start, x.end)):
            if exon.hasOverlap(curr_domain):
                segments.append("{}-{}".format(
                    max(curr_domain.start, exon.start),
                    min(curr_domain.end, exon.end)
                ))
        if len(segments) > 1:
            curr_domain.annot["sub_segments"] = ",".join(segments)

    # Write output
    log.info("Write output.")
    with GFF3IO(args.output_annotations, "w") as writer:
        for curr_domain in domains:
            writer.write(curr_domain)
    log.info("End of job")
