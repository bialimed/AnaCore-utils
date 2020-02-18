#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2020 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'

import os
import uuid
import tempfile
import unittest
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.dirname(CURRENT_DIR)
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']


########################################################################
#
# FUNCTIONS
#
########################################################################
class mergeVCFFusionsCallers(unittest.TestCase):
    def setUp(self):
        tmp_folder = tempfile.gettempdir()
        unique_id = str(uuid.uuid1())

        # Temporary files
        self.tmp_in_manta = os.path.join(tmp_folder, unique_id + "_manta.vcf")
        self.tmp_in_starfusion = os.path.join(tmp_folder, unique_id + "_starfusion.vcf")
        self.tmp_in_arriba = os.path.join(tmp_folder, unique_id + "_arriba.vcf")
        self.tmp_output = os.path.join(tmp_folder, unique_id + "_out.vcf")

        # Exec command
        self.cmd = [
            "mergeVCFFusionsCallers.py",
            "--logging-level", "DEBUG",
            "--annotation-field", "TESTANN",
            "--calling-sources", "manta", "starfusion", "arriba",
            "--inputs-variants", self.tmp_in_manta, self.tmp_in_starfusion, self.tmp_in_arriba,
            "--output-variants", self.tmp_output
        ]

        # Create input manta
        content = """##fileformat=VCFv4.1
##fileDate=20200202
##source=GenerateSVCandidates 1.6.0
##reference=file://Homo_sapiens.GRCh38.94.dna.fa
##contig=<ID=1,length=248956422>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=2,length=242193529>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=MT,length=16569>
##contig=<ID=X,length=156040895>
##contig=<ID=Y,length=57227415>
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="CIGAR alignment for each alternate indel allele">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical homology at event breakpoints">
##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description="Length of insertion">
##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">
##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length">
##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length">
##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local translocation breakend">
##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote translocation mate breakend">
##INFO=<ID=REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at this breakend">
##INFO=<ID=MATE_REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at the other breakend">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this break-end is 5' in the fusion transcript">
##INFO=<ID=RNA_STRANDED,Number=0,Type=Flag,Description="For RNA fusions, the direction of transcription is known">
##INFO=<ID=RNA_FwRvReads,Number=2,Type=Integer,Description="For RNA fusions, number of stranded reads supporting forward or reverse direction of transcription">
##INFO=<ID=RNA_Reads,Number=1,Type=Integer,Description="The number of reads and pairs that potentially support this candidate before refinement and scoring">
##INFO=<ID=RNA_CONTIG,Number=1,Type=String,Description="The sequence of the breakend spanning contig">
##INFO=<ID=RNA_CONTIG_ALN,Number=2,Type=Integer,Description="Length of the spanning contig alignment on each breakend">
##FORMAT=<ID=PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed">
##FILTER=<ID=LowEvidence,Description="RNA fusion calls without both split read and spanning pair support">
##FILTER=<ID=Imprecise,Description="RNA fusion candidates for which no spanning contig was found">
##FILTER=<ID=Local,Description="RNA call covering short genomic distance">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##cmdline=configManta.py --rna --config structural_variants/manta/splA_manta/configManta.py.ini --referenceFasta Homo_sapiens.GRCh38.94.dna.fa --bam splA_Aligned.sortedByCoord.bam --runDir splA_manta
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
1	154944376	MantaBND:20:0:1:0:0:0:1	A	]7:72948502]A	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:20:0:1:0:0:0:0;IMPRECISE;CIPOS=-81,82;BND_DEPTH=0;MATE_BND_DEPTH=78;REF_COUNT=0;MATE_REF_COUNT=6;RNA_FIRST;RNA_STRANDED	PR	6,3
17	20204330	MantaBND:120:0:1:0:0:0:0	A	]17:20205909]A	.	Local	SVTYPE=BND;MATEID=MantaBND:120:0:1:0:0:0:1;CIPOS=0,4;HOMLEN=4;HOMSEQ=AGGG;BND_DEPTH=200;MATE_BND_DEPTH=248;REF_COUNT=74;MATE_REF_COUNT=26;RNA_STRANDED;RNA_FwRvReads=18,0;RNA_Reads=23;RNA_CONTIG=AGCTCAGACAAGAATTACTAAAGGCAAACGGTGAAATTAAACATGTTTCCAGTCTGCTGGCCAAGGGGCCTTTACAACAACTAAACGGACAGGCATTCCAGCCCCACGGGAATTTTCAGTAACTGTC;RNA_CONTIG_ALN=62,65	PR:SR	108,18:14,13
17	20205909	MantaBND:120:0:1:0:0:0:1	C	C[17:20204330[	.	Local	SVTYPE=BND;MATEID=MantaBND:120:0:1:0:0:0:0;CIPOS=0,4;HOMLEN=4;HOMSEQ=AAGG;BND_DEPTH=248;MATE_BND_DEPTH=200;REF_COUNT=26;MATE_REF_COUNT=74;RNA_FIRST;RNA_STRANDED	PR:SR	108,18:14,13
2	54668542	MantaBND:196:0:2:0:0:0:0	G	]5:79693730]G	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:196:0:2:0:0:0:1;IMPRECISE;CIPOS=-82,82;BND_DEPTH=593;MATE_BND_DEPTH=0;REF_COUNT=161;MATE_REF_COUNT=0;RNA_STRANDED;RNA_FwRvReads=3,0;RNA_Reads=6;RNA_CONTIG=AAGTGGGAAAAGGACAAAGAGAAAGACAAAGAGAAGCGGTTCAGCCTTTTTGGCAAAAAGAAATGAACTCCTTTCCTTCACCTCCTGCCCTTCTCTTACCTTTTCAGTGAAATT;RNA_CONTIG_ALN=0,114	PR	162,3
20	41421600	MantaBND:215:0:1:0:0:0:1	A	]4:173113703]A	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:215:0:1:0:0:0:0;IMPRECISE;CIPOS=-80,80;BND_DEPTH=19;MATE_BND_DEPTH=0;REF_COUNT=3;MATE_REF_COUNT=0;RNA_FIRST;RNA_STRANDED	PR	3,3
22	24334571	MantaBND:244:0:1:0:0:0:1	A	A[9:84867241[	.	PASS	SVTYPE=BND;MATEID=MantaBND:244:0:1:0:0:0:0;CIPOS=0,3;HOMLEN=3;HOMSEQ=AGG;BND_DEPTH=107;MATE_BND_DEPTH=115;REF_COUNT=0;MATE_REF_COUNT=21;RNA_FIRST;RNA_STRANDED	PR:SR	7,3:15,20
4	173113703	MantaBND:215:0:1:0:0:0:0	G	G[20:41421600[	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:215:0:1:0:0:0:1;IMPRECISE;CIPOS=-80,81;BND_DEPTH=0;MATE_BND_DEPTH=19;REF_COUNT=0;MATE_REF_COUNT=3;RNA_STRANDED;RNA_FwRvReads=3,0;RNA_Reads=6	PR	3,3
5	79693730	MantaBND:196:0:2:0:0:0:1	G	G[2:54668542[	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:196:0:2:0:0:0:0;IMPRECISE;CIPOS=-81,82;BND_DEPTH=0;MATE_BND_DEPTH=593;REF_COUNT=0;MATE_REF_COUNT=161;RNA_FIRST;RNA_STRANDED	PR	162,3
6	7181123	MantaBND:313:0:1:0:0:0:0	G	]6:7189322]G	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:313:0:1:0:0:0:1;IMPRECISE;CIPOS=-19,20;BND_DEPTH=61;MATE_BND_DEPTH=62;REF_COUNT=6;MATE_REF_COUNT=10;RNA_STRANDED;RNA_FwRvReads=3,0;RNA_Reads=3;RNA_CONTIG=GAGGCCTTACAAGTGCACTGTGTGTGGCCAGTCATTTACCACCAATGGGAACATGCACAGTGTCAACGAGTACTACCA;RNA_CONTIG_ALN=60,0	PR	52,3
6	7189322	MantaBND:313:0:1:0:0:0:1	G	G[6:7181123[	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:313:0:1:0:0:0:0;IMPRECISE;CIPOS=-19,20;BND_DEPTH=62;MATE_BND_DEPTH=61;REF_COUNT=10;MATE_REF_COUNT=6;RNA_FIRST;RNA_STRANDED	PR	52,3
6	118511298	MantaBND:317:0:1:0:0:0:1	C	]6:118566316]C	.	Local	SVTYPE=BND;MATEID=MantaBND:317:0:1:0:0:0:0;CIPOS=0,3;HOMLEN=3;HOMSEQ=TGC;BND_DEPTH=74;MATE_BND_DEPTH=213;REF_COUNT=42;MATE_REF_COUNT=25;RNA_FIRST;RNA_STRANDED	PR:SR	57,13:27,9
6	118566316	MantaBND:317:0:1:0:0:0:0	T	T[6:118511298[	.	Local	SVTYPE=BND;MATEID=MantaBND:317:0:1:0:0:0:1;CIPOS=0,3;HOMLEN=3;HOMSEQ=CTG;BND_DEPTH=213;MATE_BND_DEPTH=74;REF_COUNT=25;MATE_REF_COUNT=42;RNA_STRANDED;RNA_FwRvReads=13,0;RNA_Reads=17;RNA_CONTIG=CAGTAGGAAGAGTAATCAATGATTGACTAGGCTTAAAAGATAATGTGCCACTTGAAGTTGAATGATCTGCAAACTAGCCACATAAGAATCCTCACAATTTACATAATGTCCTAACATGGCAT;RNA_CONTIG_ALN=66,56	PR:SR	57,13:27,9
7	27895220	MantaBND:323:0:1:0:0:0:0	C	]7:27991981]C	.	Local	SVTYPE=BND;MATEID=MantaBND:323:0:1:0:0:0:1;CIPOS=0,4;HOMLEN=4;HOMSEQ=TGTC;BND_DEPTH=25;MATE_BND_DEPTH=49;REF_COUNT=3;MATE_REF_COUNT=4;RNA_FIRST;RNA_STRANDED;RNA_FwRvReads=0,4;RNA_Reads=4;RNA_CONTIG=GTTGGCTGCTGTAATTCTTGTTTTTCTAAAACCCGTGGATCTGTATCTGTCGGAGTGCTGCTGCGGAATGAAGAGGAGGG;RNA_CONTIG_ALN=46,34	PR:SR	21,4:3,4
7	27991981	MantaBND:323:0:1:0:0:0:1	T	T[7:27895220[	.	Local	SVTYPE=BND;MATEID=MantaBND:323:0:1:0:0:0:0;CIPOS=0,4;HOMLEN=4;HOMSEQ=CTGT;BND_DEPTH=49;MATE_BND_DEPTH=25;REF_COUNT=4;MATE_REF_COUNT=3;RNA_STRANDED	PR:SR	21,4:3,4
7	72948502	MantaBND:20:0:1:0:0:0:0	T	T[1:154944376[	.	Imprecise	SVTYPE=BND;MATEID=MantaBND:20:0:1:0:0:0:1;IMPRECISE;CIPOS=-80,81;BND_DEPTH=78;MATE_BND_DEPTH=0;REF_COUNT=6;MATE_REF_COUNT=0;RNA_STRANDED;RNA_FwRvReads=3,0;RNA_Reads=6	PR	6,3
9	84867241	MantaBND:244:0:1:0:0:0:0	A	]22:24334571]A	.	PASS	SVTYPE=BND;MATEID=MantaBND:244:0:1:0:0:0:1;CIPOS=0,3;HOMLEN=3;HOMSEQ=GGC;BND_DEPTH=115;MATE_BND_DEPTH=107;REF_COUNT=21;MATE_REF_COUNT=0;RNA_STRANDED;RNA_FwRvReads=22,0;RNA_Reads=25;RNA_CONTIG=TCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAGGCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATCACATCTCCAAT;RNA_CONTIG_ALN=65,64	PR:SR	7,3:15,20"""
        with open(self.tmp_in_manta, "w") as writer:
            writer.write(content)

        # Create input starfusion
        content = """##fileformat=VCFv4.1
##INFO=<ID=BREAK_DINUC,Number=1,Type=String,Description="Dinucleotides flanking the breakend of the fragment excluded by the fusion.">
##INFO=<ID=BREAK_ENTROPY,Number=1,Type=Float,Description="Shannon entropy of the 15 exonic bases flanking the breakpoint. The maximum entropy is 2, representing highest complexity. The lowest would be zero (involving a 15 base mononucleotide run). Low entropy sites should generally be treated as less confident breakpoints.">
##INFO=<ID=FCANN,Number=.,Type=String,Description="Annotation generated by FusionAnnotator (see: https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki). Format: SYMBOL|Gene|Tags">
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this breakend is 5' in the fusion transcript.">
##INFO=<ID=SPLICE_TYPE,Number=1,Type=String,Description="Whether the proposed breakpoint occurs at reference exon junctions as provided by the reference transcript structure annotations (ex. gencode).">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##FORMAT=<ID=FFPM,Number=1,Type=Float,Description="Normalized measures of the fusion-supporting rna-seq fragments (fusion fragments per million total reads).">
##FORMAT=<ID=JRL,Number=1,Type=String,Description="RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (SpanningFragCount).">
##FORMAT=<ID=SFL,Number=1,Type=String,Description="RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction (JunctionReadCount).">
##FORMAT=<ID=hasLAS,Number=1,Type=String,Description="This column indicates whether there are split reads that provide 'long' (set to length of 25 bases) alignments on both sides of the putative breakpoint (LargeAnchorSupport).">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
22	24334573	2d095596-d654-47e4-bb01-f6cc093dcf4a	N	N[9:84867243[	.		BREAK_DINUC=GT;BREAK_ENTROPY=1.9656;FCANN=SPECC1L-ADORA2A|ENSG00000258555.6|INTERCHROMOSOMAL[22--9];MATEID=e83e61aa-e2bb-4aef-8245-7534dee3d62d;RNA_FIRST;SPLICE_TYPE=ONLY_REF_SPLICE;SVTYPE=BND	FFPM:JRL:PR:SFL:SR:hasLAS	8.3739:M70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A10446%3A14926%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1113%3A28221%3A16893%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A8397%3A19685%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2103%3A12190%3A5453%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1113%3A23518%3A22366%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A18236%3A20800%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A16454%3A7184%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1101%3A9462%3A22886%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2108%3A9880%3A20539%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1102%3A20172%3A13788%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1102%3A11424%3A4669%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A8409%3A3184%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A18395%3A4305%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2113%3A3241%3A11798%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1109%3A21996%3A25012%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1112%3A11576%3A18265%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A21112%3A5976%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A14979%3A6461%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2107%3A25495%3A22290%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A24296%3A7384%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A13337%3A15870%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A5888%3A7340:6:M70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A12031%3A4135%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A25374%3A7116%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A17303%3A2912%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A14548%3A7459%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A17060%3A1605%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A20858%3A19772:22:1
9	84867243	e83e61aa-e2bb-4aef-8245-7534dee3d62d	N	]22:24334573]N	.		BREAK_DINUC=AG;BREAK_ENTROPY=1.8295;FCANN=NTRK2|ENSG00000148053.16|INTERCHROMOSOMAL[22--9];MATEID=2d095596-d654-47e4-bb01-f6cc093dcf4a;SPLICE_TYPE=ONLY_REF_SPLICE;SVTYPE=BND	FFPM:PR:SR:hasLAS	8.3739:6:22:1
22	24334573	230637aa-0101-42cd-82ba-56dec3b44671	N	N[9:84861040[	.		BREAK_DINUC=GT;BREAK_ENTROPY=1.9656;FCANN=SPECC1L-ADORA2A|ENSG00000258555.6|INTERCHROMOSOMAL[22--9];MATEID=96ac47f1-83c1-4a9c-b07a-958e6c1d9d58;RNA_FIRST;SPLICE_TYPE=ONLY_REF_SPLICE;SVTYPE=BND	FFPM:JRL:PR:SFL:SR:hasLAS	2.6916:M70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A22561%3A22702%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1112%3A8112%3A18537%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A16071%3A11466:6:M70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A17303%3A2912%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A14548%3A7459%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A17060%3A1605%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A20858%3A19772%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A25374%3A7116%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A12031%3A4135:3:1
9	84861040	96ac47f1-83c1-4a9c-b07a-958e6c1d9d58	N	]22:24334573]N	.		BREAK_DINUC=AG;BREAK_ENTROPY=1.7232;FCANN=NTRK2|ENSG00000148053.16|INTERCHROMOSOMAL[22--9];MATEID=230637aa-0101-42cd-82ba-56dec3b44671;SPLICE_TYPE=ONLY_REF_SPLICE;SVTYPE=BND	FFPM:PR:SR:hasLAS	2.6916:6:3:1"""
        with open(self.tmp_in_starfusion, "w") as writer:
            writer.write(content)

        # Create input arriba
        content = """##fileformat=VCFv4.1
##INFO=<ID=GBP,Number=1,Type=String,Description="The coordinates of the genomic breakpoint which is closest to the transcriptomic breakpoint.">
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=RNA_CONTIG,Number=1,Type=String,Description="The transcript sequence assembled from the supporting reads of the most highly expressed transcript.">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this break-end is 5' in the fusion transcript.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=TESTANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|STRAND|Site|Type|GENE_SHARD|FRAMESHIFT|Protein_contig">
##FORMAT=<ID=CFD,Number=1,Type=String,Description="Each prediction is assigned one of the confidences low, medium, or high. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations.">
##FORMAT=<ID=DPS,Number=1,Type=Integer,Description="Coverage near breakpoint. The coverage is calculated as the number of fragments near the breakpoint on the side of the breakpoint that is retained in the fusion transcript. Note that the coverage calculation counts all fragments (even duplicates).">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.">
##FORMAT=<ID=RFIL,Number=.,Type=String,Description="Filters which removed one or more of the supporting reads. The number of filtered reads is given in parantheses after the name of the filter. If a filter discarded the event as a whole (all reads), the number of filtered reads is missing.">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.">
##FORMAT=<ID=SR1,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on first shard.">
##FORMAT=<ID=SR2,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on second shard.">
##FORMAT=<ID=SRL,Number=.,Type=String,Description="The names of the supporting reads.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
22	24334573	78e3738c-eb86-49d7-abd8-c5edc42798c8	N	N[9:84867243[	.	PASS	MATEID=8d6af3e4-8513-4e74-9cad-32620df27637;RNA_FIRST;SVTYPE=BND;TESTANN=SPECC1L|+|splice-site|translocation|down|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:157:3:duplicates(105):19:16:3
9	84867243	8d6af3e4-8513-4e74-9cad-32620df27637	N	]22:24334573]N	.	PASS	MATEID=78e3738c-eb86-49d7-abd8-c5edc42798c8;RNA_CONTIG=AATAAAGTCACGCAA___GCAAGAGGAGGAGCGAGGCCGGGTATACAATTACATGAATGCCGTTGAGAGAGATTTGGCAGCCTTAAGGCAGGGAATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATCACATCTCCAATGGGAGTAACACTCCATCTTCTTCGGAAGGTGGCCCAGATGCTGTCATTATTGGAATGACCAAGATCCCTGTCATTGAAAATCCC;SVTYPE=BND;TESTANN=NTRK2|+|splice-site|translocation|up|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:687:3:duplicates(105):19:16:3
22	24334573	14afe78f-c7bd-43a7-82c9-a5e0d282ee13	N	N[9:84861040[	.	PASS	MATEID=81244f10-28e7-4c28-8024-89861afecff6;RNA_FIRST;SVTYPE=BND;TESTANN=SPECC1L|+|splice-site|translocation|down|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:157:6:duplicates(127):3:0:3
9	84861040	81244f10-28e7-4c28-8024-89861afecff6	N	]22:24334573]N	.	PASS	MATEID=14afe78f-c7bd-43a7-82c9-a5e0d282ee13;RNA_CONTIG=ATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@ATTTCTCATGGTTTGGATTTGGGAAAGTAAAATCAAGACAAGGTGTTG___GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATC;SVTYPE=BND;TESTANN=NTRK2|+|splice-site|translocation|up|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:200:6:duplicates(127):3:0:3
22	24334573	a4b58e93-3461-4f31-ad51-31800be86928	N	N[9:84867243[	.	PASS	MATEID=49bfd3ed-8f18-4a4b-a26e-953454d0bf6c;RNA_FIRST;SVTYPE=BND;TESTANN=SPECC1L-ADORA2A|+|splice-site|translocation|down|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:157:3:duplicates(105):19:16:3
9	84867243	49bfd3ed-8f18-4a4b-a26e-953454d0bf6c	N	]22:24334573]N	.	PASS	MATEID=a4b58e93-3461-4f31-ad51-31800be86928;RNA_CONTIG=AATAAAGTCACGCAA___GCAAGAGGAGGAGCGAGGCCGGGTATACAATTACATGAATGCCGTTGAGAGAGATTTGGCAGCCTTAAGGCAGGGAATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATCACATCTCCAATGGGAGTAACACTCCATCTTCTTCGGAAGGTGGCCCAGATGCTGTCATTATTGGAATGACCAAGATCCCTGTCATTGAAAATCCC;SVTYPE=BND;TESTANN=NTRK2|+|splice-site|translocation|up|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:687:3:duplicates(105):19:16:3
22	24334573	8e7f19f8-43f3-49c3-89d8-4cbbb41fcbe0	N	N[9:84861040[	.	PASS	MATEID=d6ce0eef-a08e-4f59-9cf8-99d230a1863e;RNA_FIRST;SVTYPE=BND;TESTANN=SPECC1L-ADORA2A|+|splice-site|translocation|down|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:157:6:duplicates(127):3:0:3
9	84861040	d6ce0eef-a08e-4f59-9cf8-99d230a1863e	N	]22:24334573]N	.	PASS	MATEID=8e7f19f8-43f3-49c3-89d8-4cbbb41fcbe0;RNA_CONTIG=ATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@ATTTCTCATGGTTTGGATTTGGGAAAGTAAAATCAAGACAAGGTGTTG___GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATC;SVTYPE=BND;TESTANN=NTRK2|+|splice-site|translocation|up|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	CFD:DPS:PR:RFIL:SR:SR1:SR2	high:200:6:duplicates(127):3:0:3
21	8401894	2baa916b-e203-4f2a-a460-8f72c77d9c81	N	N]7:75198002]	.	PASS	MATEID=be10983b-5879-44ba-9191-3dd591491524;RNA_FIRST;SVTYPE=BND;TESTANN=FP236383.1|+|intron|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	medium:114:1:duplicates(34):3:0:3
7	75198002	be10983b-5879-44ba-9191-3dd591491524	N	N]21:8401894]	.	PASS	MATEID=2baa916b-e203-4f2a-a460-8f72c77d9c81;RNA_CONTIG=GGCGGCCGCCCCCTCGCCCGTCACGCACCGCACGTTCGTGGGGAACCTGGCGCTAAACCATTCGTAGACGACC@CAGACTAACGGTTCTAACGTTCCCTTCAAGCCACGAGGGAGAGAGTTTTCCTTTG___AGGCCTGGAATGCCAAAATCACGGAC;SVTYPE=BND;TESTANN=GTF2IP1|-|exon|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	medium:91:1:duplicates(34):3:0:3
15	33590000	4a11b081-edb8-4125-a608-67a86ec626ac	N	N[2:32441423[	.	PASS	MATEID=adbb5962-4b8c-48ee-b74c-2bb394706d6c;RNA_FIRST;SVTYPE=BND;TESTANN=RYR3|+|intron|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	medium:0:1:duplicates(18):1:0:1
2	32441423	adbb5962-4b8c-48ee-b74c-2bb394706d6c	N	]15:33590000]N	.	PASS	MATEID=4a11b081-edb8-4125-a608-67a86ec626ac;RNA_CONTIG=TGATTTCACCTGGGTGCAGGTGGGCTAAGTTCGAAAAGAGAGTCA@CCATTCGAAGATTTAAGAAAACCTCAATTTC...AGAAAAACAGATGATGGCCAGATCACAGAACATGCCCAGAGCCTTGTGTTGGATACTCTCTGTTGGTTAGCTGGAG;SVTYPE=BND;TESTANN=BIRC6|+|CDS|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	medium:298:1:duplicates(18):1:0:1
9	131198324	fbd597c8-faab-4ede-af7e-d09cfb12e5fe	N	N[21:8258375[	.	PASS	MATEID=d6d994a3-1e46-46a0-9de8-e6823ba9a62d;RNA_FIRST;SVTYPE=BND;TESTANN=NUP214|+|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:3990:7:duplicates(26):0:0:0
21	8258375	d6d994a3-1e46-46a0-9de8-e6823ba9a62d	N	]9:131198324]N	.	PASS	MATEID=fbd597c8-faab-4ede-af7e-d09cfb12e5fe;RNA_CONTIG=CTTTTCTGTGCCTGGGCAGACTGCTGTCACAGCAGCTGCTATCTCAAGTGCAGGCCCTGTGGCCGTCGAA...@...GGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGT;SVTYPE=BND;TESTANN=RNA5-8SN2(1442)%2CFP236383.1(122290)|.|intergenic|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:636:7:duplicates(26):0:0:0
1	154944457	3dad09b4-80ac-47af-b398-324f0417f97a	N	]7:72948421]N	.	PASS	MATEID=0eb6d268-a5ea-4982-8ebd-5bbfaf1f50c0;RNA_FIRST;SVTYPE=BND;TESTANN=PBXIP1|-|3'UTR|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:7:3:duplicates(4):0:0:0
7	72948421	0eb6d268-a5ea-4982-8ebd-5bbfaf1f50c0	N	N[1:154944457[	.	PASS	MATEID=3dad09b4-80ac-47af-b398-324f0417f97a;RNA_CONTIG=CCCCCCTTTTTTTTTTTACCCCTGCTTCTCCCACGGCTTCACC?CCCTATGTGAACTGTAGACTCAGATCCCAATAA...@...CCCTGGGCCCTCACCGCAGGCAGCAGTTTGCGTTTTGAAAGGTTATTGGGTCCCTTCCTCGGGCTGTGTTCTTGCT;SVTYPE=BND;TESTANN=NSUN5P2|-|exon|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:158:3:duplicates(4):0:0:0
1	154944457	c2fce270-d20d-4d02-87b9-3e0611299179	N	]7:72948421]N	.	PASS	MATEID=cccb3bbc-b1ba-4efa-83bd-297f3763b992;RNA_FIRST;SVTYPE=BND;TESTANN=PBXIP1|-|3'UTR|translocation/5'-5'|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:7:3:duplicates(4):0:0:0
7	72948421	cccb3bbc-b1ba-4efa-83bd-297f3763b992	N	N[1:154944457[	.	PASS	MATEID=c2fce270-d20d-4d02-87b9-3e0611299179;RNA_CONTIG=CCCCCCTTTTTTTTTTTACCCCTGCTTCTCCCACGGCTTCACC?CCCTATGTGAACTGTAGACTCAGATCCCAATAA...@...CCCTGGGCCCTCACCGCAGGCAGCAGTTTGCGTTTTGAAAGGTTATTGGGTCCCTTCCTCGGGCTGTGTTCTTGCT;SVTYPE=BND;TESTANN=POM121|+|CDS|translocation/5'-5'|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:158:3:duplicates(4):0:0:0
16	3728788	7acaedb4-2af9-401b-baef-a4db7705a74f	N	[21:8258375[N	.	PASS	MATEID=344898f1-b3dd-4ad4-8c9d-0c4ed551a25b;RNA_FIRST;SVTYPE=BND;TESTANN=CREBBP|-|CDS|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:191:3:duplicates(10):0:0:0
21	8258375	344898f1-b3dd-4ad4-8c9d-0c4ed551a25b	N	[16:3728788[N	.	PASS	MATEID=7acaedb4-2af9-401b-baef-a4db7705a74f;RNA_CONTIG=ATCTCACCCAGCGCTCTGCAAGACCTGCTGCGGACCCTGAAGTCGCCCAGCTCCCCTCAGCAGCAACAGCAGGTGC...@...GGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGG;SVTYPE=BND;TESTANN=RNA5-8SN2(1442)%2CFP236383.1(122290)|.|intergenic|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:636:3:duplicates(10):0:0:0
17	29089979	26d33216-5f00-4c5c-883f-9561662f98eb	N	[21:8258336[N	.	PASS	MATEID=caa20929-a87c-4bd0-b972-a9c18e88ec0d;RNA_FIRST;SVTYPE=BND;TESTANN=MYO18A|-|CDS|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:553:3:duplicates(15):0:0:0
21	8258336	caa20929-a87c-4bd0-b972-a9c18e88ec0d	N	[17:29089979[N	.	PASS	MATEID=26d33216-5f00-4c5c-883f-9561662f98eb;RNA_CONTIG=GGTGGACAAGTCCCTGGTGAGCAGGCAGGAAGCTAAGATACGGGAGCTGGAGACACGCCTGGAGTTTGAAAGGACG...@...CGGTGTGAGGCCGGTAGCGGCCCCCGGCGCGCCGGGCCCGGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCC...CCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGC;SVTYPE=BND;TESTANN=RNA5-8SN2(1403)%2CFP236383.1(122329)|.|intergenic|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:292:3:duplicates(15):0:0:0
7	74735473	908790ce-6c7d-49b2-9dae-93c8058172b4	N	N[21:8258381[	.	PASS	MATEID=b5efe0c0-12d9-4b2c-92d8-04c31775e0b9;RNA_FIRST;SVTYPE=BND;TESTANN=GTF2I|+|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:0:3:duplicates(28):0:0:0
21	8258381	b5efe0c0-12d9-4b2c-92d8-04c31775e0b9	N	]7:74735473]N	.	PASS	MATEID=908790ce-6c7d-49b2-9dae-93c8058172b4;RNA_CONTIG=AGAA___ACATGAGCTTCTGAATTCAACACGTGAAGATTTACAGCTTGATAAGCCAGCTTCAGGAG___TAAAGGAAGAATG...@...TCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGC;SVTYPE=BND;TESTANN=RNA5-8SN2(1448)%2CFP236383.1(122284)|.|intergenic|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:620:3:duplicates(28):0:0:0
6	167946889	2d656445-4786-454e-a6dc-69a659f99511	N	N[21:8258203[	.	PASS	MATEID=0d5cf2a4-4bc5-4152-a0a1-15736ebe0e94;RNA_FIRST;SVTYPE=BND;TESTANN=AFDN|+|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:831:3:duplicates(31):0:0:0
21	8258203	0d5cf2a4-4bc5-4152-a0a1-15736ebe0e94	N	]6:167946889]N	.	PASS	MATEID=2d656445-4786-454e-a6dc-69a659f99511;RNA_CONTIG=GCAGAATATAGTGAACCAAAGAAATTGCCTGGTGATGACAGACTGATGAAAAATAGAGCTGATCACCGTTCCAGCC...@...GTGAACAGGGAAGAGCCCAGCGCCGAATCCCCGCCCCGCGGCGGGGCGCGGGACATGTGGCGTACGGAAGACCCGC...TCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACG;SVTYPE=BND;TESTANN=RNA5-8SN2(1270)%2CFP236383.1(122462)|.|intergenic|translocation|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:77:3:duplicates(31):0:0:0
X	118791837	dca7f41e-6ab1-4bea-a1aa-6ac00f10e7f9	N	N]7:140777046]	.	PASS	MATEID=090ce3e8-3ac4-44fd-a076-28cbca7de98a;RNA_FIRST;SVTYPE=BND;TESTANN=IL13RA1|+|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:9:2:duplicates(7):0:0:0
7	140777046	090ce3e8-3ac4-44fd-a076-28cbca7de98a	N	N]X:118791837]	.	PASS	MATEID=dca7f41e-6ab1-4bea-a1aa-6ac00f10e7f9;RNA_CONTIG=CTGGAAGAAGTACGACATCTATGAGAAGCAAACCAAGGAGGAAACCGACTCTGTAGTGCTGATAGAAAACCTGA...@...CACAAAGCCACAACTGGCTATTGTTACCCAtTGGTGTGAGGGCTCCAGCTTGTATCACCATCTCCATATCATTGAGA;SVTYPE=BND;TESTANN=BRAF|-|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:220:2:duplicates(7):0:0:0
3	186786001	2d2ff7c0-6198-4c84-8be5-00a66e0e7488	N	N]13:67060082]	.	PASS	MATEID=85ead32d-fecf-46bd-b3be-e3e652ed54d9;RNA_FIRST;SVTYPE=BND;TESTANN=EIF4A2|+|CDS|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:3040:2:duplicates(9):0:0:0
13	67060082	85ead32d-fecf-46bd-b3be-e3e652ed54d9	N	N]3:186786001]	.	PASS	MATEID=2d2ff7c0-6198-4c84-8be5-00a66e0e7488;RNA_CONTIG=CTTGTCATGCCTGCATTGGTGGAACAAATGTTCGAAATGAAATGCAAAAACTGCAGGCTGAAGCACCACATATTGT...@...GCCCAGAGTCACACAGACAGCGATAGCACAACTAGAATTCTTATCCAGGTCTCTCTGAATTTAGAGCCACACCCC;SVTYPE=BND;TESTANN=PCDH9|-|intron|translocation|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:1:2:duplicates(9):0:0:0
4	71568735	e04656ba-544e-4b6a-91c1-58a1242f7d16	N	N[4:1806658[	.	PASS	MATEID=16fcd5c9-4c31-4a04-a78a-7544793344ad;RNA_FIRST;SVTYPE=BND;TESTANN=SLC4A4|+|3'UTR|duplication|down|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:5:2:duplicates(9):0:0:0
4	1806658	16fcd5c9-4c31-4a04-a78a-7544793344ad	N	]4:71568735]N	.	PASS	MATEID=e04656ba-544e-4b6a-91c1-58a1242f7d16;RNA_CONTIG=TTTTTTCCTGCAGCAGGAAACATAGTTTTGAGTAGTTCTACCTCTTATTTGTAGCTGCCAGGCTTTCTGTAAAAAT...@...AAGCCCGCCAACTGCACACACGACCT___GTACATGATCATGCGGGAGTGCTGGCATGCCGCGCCCTCCCAGAGGCCC;SVTYPE=BND;TESTANN=FGFR3|+|CDS|duplication|up|.|.	CFD:DPS:PR:RFIL:SR:SR1:SR2	low:825:2:duplicates(9):0:0:0"""
        with open(self.tmp_in_arriba, "w") as writer:
            writer.write(content)


    def tearDown(self):
        # Clean temporary files
        for curr_file in [self.tmp_in_arriba, self.tmp_in_manta, self.tmp_in_starfusion, self.tmp_output]:
            if os.path.exists(curr_file):
                os.remove(curr_file)


    def testResults(self):
        # Execute command
        subprocess.check_call(self.cmd, stderr=subprocess.DEVNULL)

        # Validate results
        expected = """##fileformat=VCFv4.1
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=IDSRC,Number=.,Type=String,Description="ID of breakend by source">
##INFO=<ID=MATEID,Number=A,Type=String,Description="ID of mate breakend.">
##INFO=<ID=RNA_FIRST,Number=0,Type=Flag,Description="For RNA fusions, this break-end is 5' in the fusion transcript.">
##INFO=<ID=SRC,Number=.,Type=String,Description="Fusions callers where the breakend is identified. Possible values: {'manta': 's0', 'starfusion': 's1', 'arriba': 's2'}">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant.">
##INFO=<ID=s0_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local translocation breakend",Source="manta">
##INFO=<ID=s0_CIEND,Number=2,Type=Integer,Description="Confidence interval around END",Source="manta">
##INFO=<ID=s0_CIGAR,Number=A,Type=String,Description="CIGAR alignment for each alternate indel allele",Source="manta">
##INFO=<ID=s0_CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS",Source="manta">
##INFO=<ID=s0_END,Number=1,Type=Integer,Description="End position of the variant described in this record",Source="manta">
##INFO=<ID=s0_EVENT,Number=1,Type=String,Description="ID of event associated to breakend",Source="manta">
##INFO=<ID=s0_HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints",Source="manta">
##INFO=<ID=s0_HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical homology at event breakpoints",Source="manta">
##INFO=<ID=s0_IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation",Source="manta">
##INFO=<ID=s0_LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length",Source="manta">
##INFO=<ID=s0_MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote translocation mate breakend",Source="manta">
##INFO=<ID=s0_MATE_REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at the other breakend",Source="manta">
##INFO=<ID=s0_REF_COUNT,Number=1,Type=Integer,Description="The number of reads supporting the reference allele at this breakend",Source="manta">
##INFO=<ID=s0_RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length",Source="manta">
##INFO=<ID=s0_RNA_CONTIG,Number=1,Type=String,Description="The sequence of the breakend spanning contig",Source="manta">
##INFO=<ID=s0_RNA_CONTIG_ALN,Number=2,Type=Integer,Description="Length of the spanning contig alignment on each breakend",Source="manta">
##INFO=<ID=s0_RNA_FwRvReads,Number=2,Type=Integer,Description="For RNA fusions, number of stranded reads supporting forward or reverse direction of transcription",Source="manta">
##INFO=<ID=s0_RNA_Reads,Number=1,Type=Integer,Description="The number of reads and pairs that potentially support this candidate before refinement and scoring",Source="manta">
##INFO=<ID=s0_RNA_STRANDED,Number=0,Type=Flag,Description="For RNA fusions, the direction of transcription is known",Source="manta">
##INFO=<ID=s0_SVINSLEN,Number=.,Type=Integer,Description="Length of insertion",Source="manta">
##INFO=<ID=s0_SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion",Source="manta">
##INFO=<ID=s0_SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles",Source="manta">
##INFO=<ID=s0_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="manta">
##INFO=<ID=s1_BREAK_DINUC,Number=1,Type=String,Description="Dinucleotides flanking the breakend of the fragment excluded by the fusion.",Source="starfusion">
##INFO=<ID=s1_BREAK_ENTROPY,Number=1,Type=Float,Description="Shannon entropy of the 15 exonic bases flanking the breakpoint. The maximum entropy is 2, representing highest complexity. The lowest would be zero (involving a 15 base mononucleotide run). Low entropy sites should generally be treated as less confident breakpoints.",Source="starfusion">
##INFO=<ID=s1_FCANN,Number=.,Type=String,Description="Annotation generated by FusionAnnotator (see: https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki). Format: SYMBOL|Gene|Tags",Source="starfusion">
##INFO=<ID=s1_SPLICE_TYPE,Number=1,Type=String,Description="Whether the proposed breakpoint occurs at reference exon junctions as provided by the reference transcript structure annotations (ex. gencode).",Source="starfusion">
##INFO=<ID=s1_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="starfusion">
##INFO=<ID=s2_GBP,Number=1,Type=String,Description="The coordinates of the genomic breakpoint which is closest to the transcriptomic breakpoint.",Source="arriba">
##INFO=<ID=s2_RNA_CONTIG,Number=1,Type=String,Description="The transcript sequence assembled from the supporting reads of the most highly expressed transcript.",Source="arriba">
##INFO=<ID=s2_TESTANN,Number=.,Type=String,Description="Consequence annotations. Format: SYMBOL|STRAND|Site|Type|GENE_SHARD|FRAMESHIFT|Protein_contig",Source="arriba">
##INFO=<ID=s2_VCQUAL,Number=1,Type=Float,Description="The variant quality",Source="arriba">
##FILTER=<ID=s0_Imprecise,Description="RNA fusion candidates for which no spanning contig was found",Source="manta">
##FILTER=<ID=s0_Local,Description="RNA call covering short genomic distance",Source="manta">
##FILTER=<ID=s0_LowEvidence,Description="RNA fusion calls without both split read and spanning pair support",Source="manta">
##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Count of pairs of reads supporting the fusion">
##FORMAT=<ID=PRSRC,Number=.,Type=Integer,Description="Count of pairs of reads supporting the fusion by source">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Count of reads mapping on the fusion junction">
##FORMAT=<ID=SRSRC,Number=.,Type=Integer,Description="Count of reads mapping on the fusion junction by source">
##FORMAT=<ID=s0_PR,Number=R,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed",Source="manta">
##FORMAT=<ID=s0_SR,Number=R,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed",Source="manta">
##FORMAT=<ID=s1_FFPM,Number=1,Type=Float,Description="Normalized measures of the fusion-supporting rna-seq fragments (fusion fragments per million total reads).",Source="starfusion">
##FORMAT=<ID=s1_JRL,Number=1,Type=String,Description="RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.",Source="starfusion">
##FORMAT=<ID=s1_PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment (SpanningFragCount).",Source="starfusion">
##FORMAT=<ID=s1_SFL,Number=1,Type=String,Description="RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.",Source="starfusion">
##FORMAT=<ID=s1_SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction (JunctionReadCount).",Source="starfusion">
##FORMAT=<ID=s1_hasLAS,Number=1,Type=String,Description="This column indicates whether there are split reads that provide 'long' (set to length of 25 bases) alignments on both sides of the putative breakpoint (LargeAnchorSupport).",Source="starfusion">
##FORMAT=<ID=s2_CFD,Number=1,Type=String,Description="Each prediction is assigned one of the confidences low, medium, or high. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations.",Source="arriba">
##FORMAT=<ID=s2_DPS,Number=1,Type=Integer,Description="Coverage near breakpoint. The coverage is calculated as the number of fragments near the breakpoint on the side of the breakpoint that is retained in the fusion transcript. Note that the coverage calculation counts all fragments (even duplicates).",Source="arriba">
##FORMAT=<ID=s2_PR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments that encompass the fusion junction such that one read of the pair aligns to a different gene than the other paired-end read of that fragment.",Source="arriba">
##FORMAT=<ID=s2_RFIL,Number=.,Type=String,Description="Filters which removed one or more of the supporting reads. The number of filtered reads is given in parantheses after the name of the filter. If a filter discarded the event as a whole (all reads), the number of filtered reads is missing.",Source="arriba">
##FORMAT=<ID=s2_SR,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction.",Source="arriba">
##FORMAT=<ID=s2_SR1,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on first shard.",Source="arriba">
##FORMAT=<ID=s2_SR2,Number=1,Type=Integer,Description="Number of RNA-Seq fragments containing a read that aligns as a split read at the site of the putative fusion junction with an anchor on second shard.",Source="arriba">
##FORMAT=<ID=s2_SRL,Number=.,Type=String,Description="The names of the supporting reads.",Source="arriba">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	splA
1	154944376	MantaBND:20:0:1:0:0:0:1	A	]7:72948502]A	.	s0_Imprecise	CIPOS=-81,82;IDSRC=MantaBND%3A20%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A20%3A0%3A1%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=0;s0_CIPOS=-81,82;s0_IMPRECISE;s0_MATE_BND_DEPTH=78;s0_MATE_REF_COUNT=6;s0_REF_COUNT=0;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:6,3
1	154944457	3dad09b4-80ac-47af-b398-324f0417f97a	N	]7:72948421]N	.	PASS	IDSRC=3dad09b4-80ac-47af-b398-324f0417f97a;MATEID=0eb6d268-a5ea-4982-8ebd-5bbfaf1f50c0;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=PBXIP1|-|3'UTR|translocation|up|.|.,PBXIP1|-|3'UTR|translocation/5'-5'|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:7:3:duplicates(4):0:0:0
13	67060082	85ead32d-fecf-46bd-b3be-e3e652ed54d9	N	N]3:186786001]	.	PASS	IDSRC=85ead32d-fecf-46bd-b3be-e3e652ed54d9;MATEID=2d2ff7c0-6198-4c84-8be5-00a66e0e7488;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=CTTGTCATGCCTGCATTGGTGGAACAAATGTTCGAAATGAAATGCAAAAACTGCAGGCTGAAGCACCACATATTGT...@...GCCCAGAGTCACACAGACAGCGATAGCACAACTAGAATTCTTATCCAGGTCTCTCTGAATTTAGAGCCACACCCC;s2_TESTANN=PCDH9|-|intron|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:1:2:duplicates(9):0:0:0
15	33590000	4a11b081-edb8-4125-a608-67a86ec626ac	N	N[2:32441423[	.	PASS	IDSRC=4a11b081-edb8-4125-a608-67a86ec626ac;MATEID=adbb5962-4b8c-48ee-b74c-2bb394706d6c;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=RYR3|+|intron|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	1:1:1:1:medium:0:1:duplicates(18):1:0:1
16	3728788	7acaedb4-2af9-401b-baef-a4db7705a74f	N	[21:8258375[N	.	PASS	IDSRC=7acaedb4-2af9-401b-baef-a4db7705a74f;MATEID=344898f1-b3dd-4ad4-8c9d-0c4ed551a25b;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=CREBBP|-|CDS|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:191:3:duplicates(10):0:0:0
17	20204330	MantaBND:120:0:1:0:0:0:0	A	]17:20205909]A	.	s0_Local	CIPOS=0,4;IDSRC=MantaBND%3A120%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A120%3A0%3A1%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=200;s0_CIPOS=0,4;s0_HOMLEN=4;s0_HOMSEQ=AGGG;s0_MATE_BND_DEPTH=248;s0_MATE_REF_COUNT=26;s0_REF_COUNT=74;s0_RNA_CONTIG=AGCTCAGACAAGAATTACTAAAGGCAAACGGTGAAATTAAACATGTTTCCAGTCTGCTGGCCAAGGGGCCTTTACAACAACTAAACGGACAGGCATTCCAGCCCCACGGGAATTTTCAGTAACTGTC;s0_RNA_CONTIG_ALN=62,65;s0_RNA_FwRvReads=18,0;s0_RNA_Reads=23;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	18:13:18:13:108,18:14,13
17	20205909	MantaBND:120:0:1:0:0:0:1	C	C[17:20204330[	.	s0_Local	CIPOS=0,4;IDSRC=MantaBND%3A120%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A120%3A0%3A1%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=248;s0_CIPOS=0,4;s0_HOMLEN=4;s0_HOMSEQ=AAGG;s0_MATE_BND_DEPTH=200;s0_MATE_REF_COUNT=74;s0_REF_COUNT=26;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	18:13:18:13:108,18:14,13
17	29089979	26d33216-5f00-4c5c-883f-9561662f98eb	N	[21:8258336[N	.	PASS	IDSRC=26d33216-5f00-4c5c-883f-9561662f98eb;MATEID=caa20929-a87c-4bd0-b972-a9c18e88ec0d;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=MYO18A|-|CDS|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:553:3:duplicates(15):0:0:0
2	32441423	adbb5962-4b8c-48ee-b74c-2bb394706d6c	N	]15:33590000]N	.	PASS	IDSRC=adbb5962-4b8c-48ee-b74c-2bb394706d6c;MATEID=4a11b081-edb8-4125-a608-67a86ec626ac;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=TGATTTCACCTGGGTGCAGGTGGGCTAAGTTCGAAAAGAGAGTCA@CCATTCGAAGATTTAAGAAAACCTCAATTTC...AGAAAAACAGATGATGGCCAGATCACAGAACATGCCCAGAGCCTTGTGTTGGATACTCTCTGTTGGTTAGCTGGAG;s2_TESTANN=BIRC6|+|CDS|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	1:1:1:1:medium:298:1:duplicates(18):1:0:1
2	54668542	MantaBND:196:0:2:0:0:0:0	G	]5:79693730]G	.	s0_Imprecise	CIPOS=-82,82;IDSRC=MantaBND%3A196%3A0%3A2%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A196%3A0%3A2%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=593;s0_CIPOS=-82,82;s0_IMPRECISE;s0_MATE_BND_DEPTH=0;s0_MATE_REF_COUNT=0;s0_REF_COUNT=161;s0_RNA_CONTIG=AAGTGGGAAAAGGACAAAGAGAAAGACAAAGAGAAGCGGTTCAGCCTTTTTGGCAAAAAGAAATGAACTCCTTTCCTTCACCTCCTGCCCTTCTCTTACCTTTTCAGTGAAATT;s0_RNA_CONTIG_ALN=0,114;s0_RNA_FwRvReads=3,0;s0_RNA_Reads=6;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:162,3
20	41421600	MantaBND:215:0:1:0:0:0:1	A	]4:173113703]A	.	s0_Imprecise	CIPOS=-80,80;IDSRC=MantaBND%3A215%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A215%3A0%3A1%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=19;s0_CIPOS=-80,80;s0_IMPRECISE;s0_MATE_BND_DEPTH=0;s0_MATE_REF_COUNT=0;s0_REF_COUNT=3;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:3,3
21	8258203	0d5cf2a4-4bc5-4152-a0a1-15736ebe0e94	N	]6:167946889]N	.	PASS	IDSRC=0d5cf2a4-4bc5-4152-a0a1-15736ebe0e94;MATEID=2d656445-4786-454e-a6dc-69a659f99511;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=GCAGAATATAGTGAACCAAAGAAATTGCCTGGTGATGACAGACTGATGAAAAATAGAGCTGATCACCGTTCCAGCC...@...GTGAACAGGGAAGAGCCCAGCGCCGAATCCCCGCCCCGCGGCGGGGCGCGGGACATGTGGCGTACGGAAGACCCGC...TCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACG;s2_TESTANN=RNA5-8SN2(1270)%2CFP236383.1(122462)|.|intergenic|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:77:3:duplicates(31):0:0:0
21	8258336	caa20929-a87c-4bd0-b972-a9c18e88ec0d	N	[17:29089979[N	.	PASS	IDSRC=caa20929-a87c-4bd0-b972-a9c18e88ec0d;MATEID=26d33216-5f00-4c5c-883f-9561662f98eb;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=GGTGGACAAGTCCCTGGTGAGCAGGCAGGAAGCTAAGATACGGGAGCTGGAGACACGCCTGGAGTTTGAAAGGACG...@...CGGTGTGAGGCCGGTAGCGGCCCCCGGCGCGCCGGGCCCGGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCC...CCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGC;s2_TESTANN=RNA5-8SN2(1403)%2CFP236383.1(122329)|.|intergenic|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:292:3:duplicates(15):0:0:0
21	8258375	344898f1-b3dd-4ad4-8c9d-0c4ed551a25b	N	[16:3728788[N	.	PASS	IDSRC=344898f1-b3dd-4ad4-8c9d-0c4ed551a25b;MATEID=7acaedb4-2af9-401b-baef-a4db7705a74f;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=ATCTCACCCAGCGCTCTGCAAGACCTGCTGCGGACCCTGAAGTCGCCCAGCTCCCCTCAGCAGCAACAGCAGGTGC...@...GGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGG;s2_TESTANN=RNA5-8SN2(1442)%2CFP236383.1(122290)|.|intergenic|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:636:3:duplicates(10):0:0:0
21	8258375	d6d994a3-1e46-46a0-9de8-e6823ba9a62d	N	]9:131198324]N	.	PASS	IDSRC=d6d994a3-1e46-46a0-9de8-e6823ba9a62d;MATEID=fbd597c8-faab-4ede-af7e-d09cfb12e5fe;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=CTTTTCTGTGCCTGGGCAGACTGCTGTCACAGCAGCTGCTATCTCAAGTGCAGGCCCTGTGGCCGTCGAA...@...GGGTCTTCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGT;s2_TESTANN=RNA5-8SN2(1442)%2CFP236383.1(122290)|.|intergenic|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	7:0:7:0:low:636:7:duplicates(26):0:0:0
21	8258381	b5efe0c0-12d9-4b2c-92d8-04c31775e0b9	N	]7:74735473]N	.	PASS	IDSRC=b5efe0c0-12d9-4b2c-92d8-04c31775e0b9;MATEID=908790ce-6c7d-49b2-9dae-93c8058172b4;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=AGAA___ACATGAGCTTCTGAATTCAACACGTGAAGATTTACAGCTTGATAAGCCAGCTTCAGGAG___TAAAGGAAGAATG...@...TCCCGGAGTCGGGTTGCTTGGGAATGCAGCCCAAAGCGGGTGGTAAACTCCATCTAAGGCTAAATACCGGCACGAGACCGATAGTCAACAAGTACCGTAAGGGAAAGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGC;s2_TESTANN=RNA5-8SN2(1448)%2CFP236383.1(122284)|.|intergenic|translocation|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:620:3:duplicates(28):0:0:0
21	8401894	2baa916b-e203-4f2a-a460-8f72c77d9c81	N	N]7:75198002]	.	PASS	IDSRC=2baa916b-e203-4f2a-a460-8f72c77d9c81;MATEID=be10983b-5879-44ba-9191-3dd591491524;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=FP236383.1|+|intron|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	1:3:1:3:medium:114:1:duplicates(34):3:0:3
22	24334571	MantaBND:244:0:1:0:0:0:1	A	A[9:84867241[	.	PASS	CIPOS=0,3;IDSRC=MantaBND%3A244%3A0%3A1%3A0%3A0%3A0%3A1,2d095596-d654-47e4-bb01-f6cc093dcf4a,78e3738c-eb86-49d7-abd8-c5edc42798c8;MATEID=8d6af3e4-8513-4e74-9cad-32620df27637;RNA_FIRST;SRC=manta,starfusion,arriba;SVTYPE=BND;s0_BND_DEPTH=107;s0_CIPOS=0,3;s0_HOMLEN=3;s0_HOMSEQ=AGG;s0_MATE_BND_DEPTH=115;s0_MATE_REF_COUNT=21;s0_REF_COUNT=0;s0_RNA_STRANDED;s1_BREAK_DINUC=GT;s1_BREAK_ENTROPY=1.9656;s1_FCANN=SPECC1L-ADORA2A|ENSG00000258555.6|INTERCHROMOSOMAL[22--9];s1_SPLICE_TYPE=ONLY_REF_SPLICE;s2_TESTANN=SPECC1L|+|splice-site|translocation|down|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN,SPECC1L-ADORA2A|+|splice-site|translocation|down|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR:s1_FFPM:s1_JRL:s1_PR:s1_SFL:s1_SR:s1_hasLAS:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:20:3,6,3:20,22,19:7,3:15,20:8.3739:M70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A10446%3A14926%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1113%3A28221%3A16893%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A8397%3A19685%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2103%3A12190%3A5453%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1113%3A23518%3A22366%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A18236%3A20800%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A16454%3A7184%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1101%3A9462%3A22886%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2108%3A9880%3A20539%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1102%3A20172%3A13788%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1102%3A11424%3A4669%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A8409%3A3184%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A18395%3A4305%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2113%3A3241%3A11798%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1109%3A21996%3A25012%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1112%3A11576%3A18265%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A21112%3A5976%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A14979%3A6461%2CM70265%3A74%3A000000000-B4F4J%3A1%3A2107%3A25495%3A22290%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A24296%3A7384%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A13337%3A15870%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A5888%3A7340:6:M70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A12031%3A4135%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A25374%3A7116%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A17303%3A2912%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A14548%3A7459%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A17060%3A1605%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A20858%3A19772:22:1:high:157:3:duplicates(105):19:16:3
22	24334573	230637aa-0101-42cd-82ba-56dec3b44671	N	N[9:84861040[	.	PASS	IDSRC=230637aa-0101-42cd-82ba-56dec3b44671,14afe78f-c7bd-43a7-82c9-a5e0d282ee13;MATEID=81244f10-28e7-4c28-8024-89861afecff6;RNA_FIRST;SRC=starfusion,arriba;SVTYPE=BND;s1_BREAK_DINUC=GT;s1_BREAK_ENTROPY=1.9656;s1_FCANN=SPECC1L-ADORA2A|ENSG00000258555.6|INTERCHROMOSOMAL[22--9];s1_SPLICE_TYPE=ONLY_REF_SPLICE;s2_TESTANN=SPECC1L|+|splice-site|translocation|down|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH,SPECC1L-ADORA2A|+|splice-site|translocation|down|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	PR:SR:PRSRC:SRSRC:s1_FFPM:s1_JRL:s1_PR:s1_SFL:s1_SR:s1_hasLAS:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	6:3:6,6:3,3:2.6916:M70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A22561%3A22702%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1112%3A8112%3A18537%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A16071%3A11466:6:M70265%3A74%3A000000000-B4F4J%3A1%3A1107%3A17303%3A2912%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A14548%3A7459%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1114%3A17060%3A1605%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1106%3A20858%3A19772%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1108%3A25374%3A7116%2CM70265%3A74%3A000000000-B4F4J%3A1%3A1119%3A12031%3A4135:3:1:high:157:6:duplicates(127):3:0:3
3	186786001	2d2ff7c0-6198-4c84-8be5-00a66e0e7488	N	N]13:67060082]	.	PASS	IDSRC=2d2ff7c0-6198-4c84-8be5-00a66e0e7488;MATEID=85ead32d-fecf-46bd-b3be-e3e652ed54d9;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=EIF4A2|+|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:3040:2:duplicates(9):0:0:0
4	1806658	16fcd5c9-4c31-4a04-a78a-7544793344ad	N	]4:71568735]N	.	PASS	IDSRC=16fcd5c9-4c31-4a04-a78a-7544793344ad;MATEID=e04656ba-544e-4b6a-91c1-58a1242f7d16;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=TTTTTTCCTGCAGCAGGAAACATAGTTTTGAGTAGTTCTACCTCTTATTTGTAGCTGCCAGGCTTTCTGTAAAAAT...@...AAGCCCGCCAACTGCACACACGACCT___GTACATGATCATGCGGGAGTGCTGGCATGCCGCGCCCTCCCAGAGGCCC;s2_TESTANN=FGFR3|+|CDS|duplication|up|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:825:2:duplicates(9):0:0:0
4	71568735	e04656ba-544e-4b6a-91c1-58a1242f7d16	N	N[4:1806658[	.	PASS	IDSRC=e04656ba-544e-4b6a-91c1-58a1242f7d16;MATEID=16fcd5c9-4c31-4a04-a78a-7544793344ad;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=SLC4A4|+|3'UTR|duplication|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:5:2:duplicates(9):0:0:0
4	173113703	MantaBND:215:0:1:0:0:0:0	G	G[20:41421600[	.	s0_Imprecise	CIPOS=-80,81;IDSRC=MantaBND%3A215%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A215%3A0%3A1%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=0;s0_CIPOS=-80,81;s0_IMPRECISE;s0_MATE_BND_DEPTH=19;s0_MATE_REF_COUNT=3;s0_REF_COUNT=0;s0_RNA_FwRvReads=3,0;s0_RNA_Reads=6;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:3,3
5	79693730	MantaBND:196:0:2:0:0:0:1	G	G[2:54668542[	.	s0_Imprecise	CIPOS=-81,82;IDSRC=MantaBND%3A196%3A0%3A2%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A196%3A0%3A2%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=0;s0_CIPOS=-81,82;s0_IMPRECISE;s0_MATE_BND_DEPTH=593;s0_MATE_REF_COUNT=161;s0_REF_COUNT=0;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:162,3
6	7181123	MantaBND:313:0:1:0:0:0:0	G	]6:7189322]G	.	s0_Imprecise	CIPOS=-19,20;IDSRC=MantaBND%3A313%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A313%3A0%3A1%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=61;s0_CIPOS=-19,20;s0_IMPRECISE;s0_MATE_BND_DEPTH=62;s0_MATE_REF_COUNT=10;s0_REF_COUNT=6;s0_RNA_CONTIG=GAGGCCTTACAAGTGCACTGTGTGTGGCCAGTCATTTACCACCAATGGGAACATGCACAGTGTCAACGAGTACTACCA;s0_RNA_CONTIG_ALN=60,0;s0_RNA_FwRvReads=3,0;s0_RNA_Reads=3;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:52,3
6	7189322	MantaBND:313:0:1:0:0:0:1	G	G[6:7181123[	.	s0_Imprecise	CIPOS=-19,20;IDSRC=MantaBND%3A313%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A313%3A0%3A1%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=62;s0_CIPOS=-19,20;s0_IMPRECISE;s0_MATE_BND_DEPTH=61;s0_MATE_REF_COUNT=6;s0_REF_COUNT=10;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:52,3
6	118511298	MantaBND:317:0:1:0:0:0:1	C	]6:118566316]C	.	s0_Local	CIPOS=0,3;IDSRC=MantaBND%3A317%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A317%3A0%3A1%3A0%3A0%3A0%3A0;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=74;s0_CIPOS=0,3;s0_HOMLEN=3;s0_HOMSEQ=TGC;s0_MATE_BND_DEPTH=213;s0_MATE_REF_COUNT=25;s0_REF_COUNT=42;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	13:9:13:9:57,13:27,9
6	118566316	MantaBND:317:0:1:0:0:0:0	T	T[6:118511298[	.	s0_Local	CIPOS=0,3;IDSRC=MantaBND%3A317%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A317%3A0%3A1%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=213;s0_CIPOS=0,3;s0_HOMLEN=3;s0_HOMSEQ=CTG;s0_MATE_BND_DEPTH=74;s0_MATE_REF_COUNT=42;s0_REF_COUNT=25;s0_RNA_CONTIG=CAGTAGGAAGAGTAATCAATGATTGACTAGGCTTAAAAGATAATGTGCCACTTGAAGTTGAATGATCTGCAAACTAGCCACATAAGAATCCTCACAATTTACATAATGTCCTAACATGGCAT;s0_RNA_CONTIG_ALN=66,56;s0_RNA_FwRvReads=13,0;s0_RNA_Reads=17;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	13:9:13:9:57,13:27,9
6	167946889	2d656445-4786-454e-a6dc-69a659f99511	N	N[21:8258203[	.	PASS	IDSRC=2d656445-4786-454e-a6dc-69a659f99511;MATEID=0d5cf2a4-4bc5-4152-a0a1-15736ebe0e94;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=AFDN|+|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:831:3:duplicates(31):0:0:0
7	27895220	MantaBND:323:0:1:0:0:0:0	C	]7:27991981]C	.	s0_Local	CIPOS=0,4;IDSRC=MantaBND%3A323%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A323%3A0%3A1%3A0%3A0%3A0%3A1;RNA_FIRST;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=25;s0_CIPOS=0,4;s0_HOMLEN=4;s0_HOMSEQ=TGTC;s0_MATE_BND_DEPTH=49;s0_MATE_REF_COUNT=4;s0_REF_COUNT=3;s0_RNA_CONTIG=GTTGGCTGCTGTAATTCTTGTTTTTCTAAAACCCGTGGATCTGTATCTGTCGGAGTGCTGCTGCGGAATGAAGAGGAGGG;s0_RNA_CONTIG_ALN=46,34;s0_RNA_FwRvReads=0,4;s0_RNA_Reads=4;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	4:4:4:4:21,4:3,4
7	27991981	MantaBND:323:0:1:0:0:0:1	T	T[7:27895220[	.	s0_Local	CIPOS=0,4;IDSRC=MantaBND%3A323%3A0%3A1%3A0%3A0%3A0%3A1;MATEID=MantaBND%3A323%3A0%3A1%3A0%3A0%3A0%3A0;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=49;s0_CIPOS=0,4;s0_HOMLEN=4;s0_HOMSEQ=CTGT;s0_MATE_BND_DEPTH=25;s0_MATE_REF_COUNT=3;s0_REF_COUNT=4;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR	4:4:4:4:21,4:3,4
7	72948421	0eb6d268-a5ea-4982-8ebd-5bbfaf1f50c0	N	N[1:154944457[	.	PASS	IDSRC=0eb6d268-a5ea-4982-8ebd-5bbfaf1f50c0;MATEID=3dad09b4-80ac-47af-b398-324f0417f97a;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=CCCCCCTTTTTTTTTTTACCCCTGCTTCTCCCACGGCTTCACC?CCCTATGTGAACTGTAGACTCAGATCCCAATAA...@...CCCTGGGCCCTCACCGCAGGCAGCAGTTTGCGTTTTGAAAGGTTATTGGGTCCCTTCCTCGGGCTGTGTTCTTGCT;s2_TESTANN=NSUN5P2|-|exon|translocation|down|.|.,POM121|+|CDS|translocation/5'-5'|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:158:3:duplicates(4):0:0:0
7	72948502	MantaBND:20:0:1:0:0:0:0	T	T[1:154944376[	.	s0_Imprecise	CIPOS=-80,81;IDSRC=MantaBND%3A20%3A0%3A1%3A0%3A0%3A0%3A0;MATEID=MantaBND%3A20%3A0%3A1%3A0%3A0%3A0%3A1;SRC=manta;SVTYPE=BND;s0_BND_DEPTH=78;s0_CIPOS=-80,81;s0_IMPRECISE;s0_MATE_BND_DEPTH=0;s0_MATE_REF_COUNT=0;s0_REF_COUNT=6;s0_RNA_FwRvReads=3,0;s0_RNA_Reads=6;s0_RNA_STRANDED	PR:SR:PRSRC:SRSRC:s0_PR	3:0:3:0:6,3
7	74735473	908790ce-6c7d-49b2-9dae-93c8058172b4	N	N[21:8258381[	.	PASS	IDSRC=908790ce-6c7d-49b2-9dae-93c8058172b4;MATEID=b5efe0c0-12d9-4b2c-92d8-04c31775e0b9;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=GTF2I|+|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:0:3:0:low:0:3:duplicates(28):0:0:0
7	75198002	be10983b-5879-44ba-9191-3dd591491524	N	N]21:8401894]	.	PASS	IDSRC=be10983b-5879-44ba-9191-3dd591491524;MATEID=2baa916b-e203-4f2a-a460-8f72c77d9c81;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=GGCGGCCGCCCCCTCGCCCGTCACGCACCGCACGTTCGTGGGGAACCTGGCGCTAAACCATTCGTAGACGACC@CAGACTAACGGTTCTAACGTTCCCTTCAAGCCACGAGGGAGAGAGTTTTCCTTTG___AGGCCTGGAATGCCAAAATCACGGAC;s2_TESTANN=GTF2IP1|-|exon|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	1:3:1:3:medium:91:1:duplicates(34):3:0:3
7	140777046	090ce3e8-3ac4-44fd-a076-28cbca7de98a	N	N]X:118791837]	.	PASS	IDSRC=090ce3e8-3ac4-44fd-a076-28cbca7de98a;MATEID=dca7f41e-6ab1-4bea-a1aa-6ac00f10e7f9;SRC=arriba;SVTYPE=BND;s2_RNA_CONTIG=CTGGAAGAAGTACGACATCTATGAGAAGCAAACCAAGGAGGAAACCGACTCTGTAGTGCTGATAGAAAACCTGA...@...CACAAAGCCACAACTGGCTATTGTTACCCAtTGGTGTGAGGGCTCCAGCTTGTATCACCATCTCCATATCATTGAGA;s2_TESTANN=BRAF|-|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:220:2:duplicates(7):0:0:0
9	84861040	96ac47f1-83c1-4a9c-b07a-958e6c1d9d58	N	]22:24334573]N	.	PASS	IDSRC=96ac47f1-83c1-4a9c-b07a-958e6c1d9d58,81244f10-28e7-4c28-8024-89861afecff6;MATEID=14afe78f-c7bd-43a7-82c9-a5e0d282ee13;SRC=starfusion,arriba;SVTYPE=BND;s1_BREAK_DINUC=AG;s1_BREAK_ENTROPY=1.7232;s1_FCANN=NTRK2|ENSG00000148053.16|INTERCHROMOSOMAL[22--9];s1_SPLICE_TYPE=ONLY_REF_SPLICE;s2_RNA_CONTIG=ATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@ATTTCTCATGGTTTGGATTTGGGAAAGTAAAATCAAGACAAGGTGTTG___GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATC;s2_TESTANN=NTRK2|+|splice-site|translocation|up|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH,NTRK2|+|splice-site|translocation|up|false|MGLSRRSSTSSEPTPTVKTLIKSFDSASQ@dFSWFGFGKVKSRQGVGPASVISNDDDSASPLH	PR:SR:PRSRC:SRSRC:s1_FFPM:s1_PR:s1_SR:s1_hasLAS:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	6:3:6,6:3,3:2.6916:6:3:1:high:200:6:duplicates(127):3:0:3
9	84867241	MantaBND:244:0:1:0:0:0:0	A	]22:24334571]A	.	PASS	CIPOS=0,3;IDSRC=MantaBND%3A244%3A0%3A1%3A0%3A0%3A0%3A0,e83e61aa-e2bb-4aef-8245-7534dee3d62d,8d6af3e4-8513-4e74-9cad-32620df27637;MATEID=78e3738c-eb86-49d7-abd8-c5edc42798c8;SRC=manta,starfusion,arriba;SVTYPE=BND;s0_BND_DEPTH=115;s0_CIPOS=0,3;s0_HOMLEN=3;s0_HOMSEQ=GGC;s0_MATE_BND_DEPTH=107;s0_MATE_REF_COUNT=0;s0_REF_COUNT=21;s0_RNA_CONTIG=TCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAGGCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATCACATCTCCAAT;s0_RNA_CONTIG_ALN=65,64;s0_RNA_FwRvReads=22,0;s0_RNA_Reads=25;s0_RNA_STRANDED;s1_BREAK_DINUC=AG;s1_BREAK_ENTROPY=1.8295;s1_FCANN=NTRK2|ENSG00000148053.16|INTERCHROMOSOMAL[22--9];s1_SPLICE_TYPE=ONLY_REF_SPLICE;s2_RNA_CONTIG=AATAAAGTCACGCAA___GCAAGAGGAGGAGCGAGGCCGGGTATACAATTACATGAATGCCGTTGAGAGAGATTTGGCAGCCTTAAGGCAGGGAATGGGACTGAGTAGAAGGTCCTCGACTTCCTCAGAGCCAACTCCTACAGTAAAAACCCTCATCAAGTCCTTTGACAGTGCATCTCAAG@GCCCAGCCTCCGTTATCAGCAATGATGATGACTCTGCCAGCCCACTCCATCACATCTCCAATGGGAGTAACACTCCATCTTCTTCGGAAGGTGGCCCAGATGCTGTCATTATTGGAATGACCAAGATCCCTGTCATTGAAAATCCC;s2_TESTANN=NTRK2|+|splice-site|translocation|up|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN,NTRK2|+|splice-site|translocation|up|false|IKSRKQEEERGRVYNYMNAVERDLAALRQGMGLSRRSSTSSEPTPTVKTLIKSFDSASQ@gPASVISNDDDSASPLHHISNGSNTPSSSEGGPDAVIIGMTKIPVIEN	PR:SR:PRSRC:SRSRC:s0_PR:s0_SR:s1_FFPM:s1_PR:s1_SR:s1_hasLAS:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	3:20:3,6,3:20,22,19:7,3:15,20:8.3739:6:22:1:high:687:3:duplicates(105):19:16:3
9	131198324	fbd597c8-faab-4ede-af7e-d09cfb12e5fe	N	N[21:8258375[	.	PASS	IDSRC=fbd597c8-faab-4ede-af7e-d09cfb12e5fe;MATEID=d6d994a3-1e46-46a0-9de8-e6823ba9a62d;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=NUP214|+|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	7:0:7:0:low:3990:7:duplicates(26):0:0:0
X	118791837	dca7f41e-6ab1-4bea-a1aa-6ac00f10e7f9	N	N]7:140777046]	.	PASS	IDSRC=dca7f41e-6ab1-4bea-a1aa-6ac00f10e7f9;MATEID=090ce3e8-3ac4-44fd-a076-28cbca7de98a;RNA_FIRST;SRC=arriba;SVTYPE=BND;s2_TESTANN=IL13RA1|+|CDS|translocation|down|.|.	PR:SR:PRSRC:SRSRC:s2_CFD:s2_DPS:s2_PR:s2_RFIL:s2_SR:s2_SR1:s2_SR2	2:0:2:0:low:9:2:duplicates(7):0:0:0"""
        observed = None
        with open(self.tmp_output) as FH_results:
            observed = FH_results.read().strip()
        # Compare
        self.assertEqual(
            sorted(expected),
            sorted(observed)
        )


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    unittest.main()
