#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2019 IUCT-O'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'escudie.frederic@iuct-oncopole.fr'
__status__ = 'prod'


from anacore.gtf import GTFIO
from anacore.sv import HashedSVIO
import argparse
import logging
import os
import re
import sys


########################################################################
#
# GENES ALIASES
#
########################################################################
def aliasesBySymbols(in_aliases):
    """
    Return all names aliases by each gene symbol.

    :param in_aliases: Path to genes synonyms (format: TSV).
    :type in_aliases: str
    :return: Names aliases by each gene symbol.
    :rtype: dict
    """
    is_ncbi = True
    with HashedSVIO(in_aliases) as reader:
        if "Gene name" in reader.titles and "Gene Synonym" in reader.titles:
            is_ncbi = False
    if is_ncbi:
        return aliasesBySymbolsFromNCBI(in_aliases)
    else:
        return aliasesBySymbolsFromEnsembl(in_aliases)


def aliasesBySymbolsFromEnsembl(in_aliases):
    """
    Return all names aliases by each gene symbol from Ensembl biomart export.

    :param in_aliases: Path to genes synonyms from Ensembl (format: TSV).
    :type in_aliases: str
    :return: Names aliases by each gene symbol.
    :rtype: dict
    """
    aliases_by_symbol = {}
    with HashedSVIO(in_aliases) as reader:
        for record in reader:
            name = record["Gene name"]
            alias = record["Gene Synonym"]
            if name not in aliases_by_symbol:
                aliases_by_symbol[name] = [name, alias]
            else:
                aliases_by_symbol[name].append(alias)
            if alias not in aliases_by_symbol:
                aliases_by_symbol[alias] = [alias, name]
            else:
                aliases_by_symbol[alias].append(name)
    return aliases_by_symbol


def aliasesBySymbolsFromNCBI(in_aliases):
    """
    Return all names aliases by each gene symbol from NCBI RefSeq gene_info.

    :param in_aliases: Path to genes synonyms from gene_info (format: TSV).
    :type in_aliases: str
    :return: Names aliases by each gene symbol.
    :rtype: dict
    """
    aliases_by_symbol = {}
    with HashedSVIO(in_aliases) as reader:
        for record in reader:
            name = record["Symbol"]
            aliases = record["Synonyms"].split("|")
            if name not in aliases_by_symbol:
                aliases_by_symbol[name] = [name] + aliases
            else:
                aliases_by_symbol[name] += aliases
            for alias in aliases:
                if alias not in aliases_by_symbol:
                    aliases_by_symbol[alias] = [name] + aliases
                else:
                    aliases_by_symbol[alias] += [name]
    return aliases_by_symbol


def annotSymbols(in_annotations):
    """
    Return list of genes names used in genes annotations file.

    :param in_annotations: Path to the genes annotations file.
    :type in_annotations: str
    :return: List of genes names used in genes annotations file.
    :rtype: set
    """
    annotation_symbols = set()
    with GTFIO(in_annotations) as reader:
        for record in reader:
            if "gene_name" in record.annot:
                annotation_symbols.add(record.annot["gene_name"])
            elif "gene" in record.annot:
                annotation_symbols.add(record.annot["gene"])
    return annotation_symbols


def selectAnnotSymbol(gene_symbol, annotation_symbols, aliases_by_symbol):
    """
    Return alias of the gene symbol used in genes annotations file.

    :param gene_symbol: Gene name.
    :type gene_symbol: str
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :return: Alias of the gene symbol used in genes annotations file.
    :rtype: str
    """
    retained_name = None
    if gene_symbol not in aliases_by_symbol and "ORF" in gene_symbol:
        gene_symbol = gene_symbol.replace("ORF", "orf")
    aliases = [gene_symbol]
    if gene_symbol in aliases_by_symbol:
        aliases = aliases_by_symbol[gene_symbol]
    for curr_name in aliases:
        if curr_name in annotation_symbols:
            retained_name = curr_name
    if retained_name is None:
        raise Exception("The gene with aliases {} cannot be found in genes annotations.".format(aliases))
    return retained_name


########################################################################
#
# FUNCTIONS
#
########################################################################
def loadBabiceanu(db_path, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols):
    """
    Set fusions partners data from single source database.

    :param db_path: Path to Babiceanu 2016 healthy fusions database (format: TSV).
    :type db_path: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (Babiceanu, BodyMap, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    """
    loadGeneric(db_path, "Babiceanu", db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols)


def loadBodyMap(db_path, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols):
    """
    Set fusions partners data from single source database.

    :param db_path: Path to BodyMap healthy fusions database (format: TSV).
    :type db_path: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (Babiceanu, BodyMap, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    """
    loadGeneric(db_path, "BodyMap", db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols)


def loadChimerdb(db_path, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols):
    """
    Set fusions partners data from chimerdb database.

    :param db_path: Path to the chimerdb database (format: TSV).
    :type db_path: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (chimerakb, cosmic, mitelman, pubmed, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    """
    # id    Source    webSource    Fusion_pair    H_gene    H_chr    H_position    H_strand    T_gene    T_chr    T_position    T_strand    Breakpoint_Type    Genome_Build_Version    PMID    Disease    Validation    Kinase    Oncogene    Tumor_suppressor    Receptor    Transcription_Factor    ChimerPub    ChimerSeq
    with HashedSVIO(db_path) as reader:
        for record in reader:
            try:
                up_gene = selectAnnotSymbol(record["H_gene"], annotation_symbols, aliases_by_symbol)
                down_gene = selectAnnotSymbol(record["T_gene"], annotation_symbols, aliases_by_symbol)
                fusion_partners = "{}_@_{}".format(up_gene, down_gene)
                source = "chimerdb_{}".format(db_version)
                if fusion_partners not in fusions_by_partners:
                    fusions_by_partners[fusion_partners] = {source: set()}
                if source not in fusions_by_partners[fusion_partners]:
                    fusions_by_partners[fusion_partners][source] = set()
                fusions_by_partners[fusion_partners][source].add(int(record["id"]))
                if record["PMID"] != "":
                    if "PMID" not in fusions_by_partners[fusion_partners]:
                        fusions_by_partners[fusion_partners]["PMID"] = set()
                    fusions_by_partners[fusion_partners]["PMID"].add(int(record["PMID"]))
            except Exception:
                print("warning", "chimeradb", db_version, record["H_gene"], record["T_gene"], sep="\t")


def loadCosmic(db_path, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols):
    """
    Set fusions partners data from cosmic database.

    :param db_path: Path to the cosmic database (format: TSV).
    :type db_path: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (chimerakb, cosmic, mitelman, pubmed, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    """
    # Sample ID    Sample name    Primary site    Site subtype 1    Site subtype 2    Site subtype 3    Primary histology    Histology subtype 1    Histology subtype 2    Histology subtype 3    Fusion ID    Translocation Name    5'_CHROMOSOME    5'_GENOME_START_FROM    5'_GENOME_START_TO    5'_GENOME_STOP_FROM    5'_GENOME_STOP_TO    5'_STRAND    3'_CHROMOSOME    3'_GENOME_START_FROM    3'_GENOME_START_TO    3'_GENOME_STOP_FROM    3'_GENOME_STOP_TO    3'_STRAND    Fusion type    Pubmed_PMID
    with HashedSVIO(db_path) as reader:
        for record in reader:
            if record["Translocation Name"] != "":
                matches = re.fullmatch(r"ENS.+\((.+)\):.+_ENS.+\((.+)\):.+", record["Translocation Name"])  # ENST00000324093.4(PLXND1):r.1_2864_ENST00000393238.3(TMCC1):r.918_5992
                if matches is None:
                    print("warning", "cosmic", db_version, record["Translocation Name"], sep="\t")
                else:
                    up_gene, down_gene = matches.groups()
                    up_gene = selectAnnotSymbol(up_gene, annotation_symbols, aliases_by_symbol)
                    down_gene = selectAnnotSymbol(down_gene, annotation_symbols, aliases_by_symbol)
                    fusion_partners = "{}_@_{}".format(up_gene, down_gene)
                    source = "cosmic_{}".format(db_version)
                    if fusion_partners not in fusions_by_partners:
                        fusions_by_partners[fusion_partners] = {source: set()}
                    if source not in fusions_by_partners[fusion_partners]:
                        fusions_by_partners[fusion_partners][source] = set()
                    fusions_by_partners[fusion_partners][source].add(int(record["Fusion ID"]))
                    if record["Pubmed_PMID"] != "":
                        if "PMID" not in fusions_by_partners[fusion_partners]:
                            fusions_by_partners[fusion_partners]["PMID"] = set()
                        fusions_by_partners[fusion_partners]["PMID"].add(int(record["Pubmed_PMID"]))


def loadGeneric(db_path, db_name, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols, up_title="up_gene", down_title="down_gene"):
    """
    Set fusions partners data from single source database.

    :param db_path: Path to fusions database (format: TSV).
    :type db_path: str
    :param db_name: Database name.
    :type db_name: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (Babiceanu, BodyMap, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    :param up_title: Title of column containing gene name of first partner.
    :type up_title: str
    :param down_title: Title of column containing gene name of second partner.
    :type down_title: str
    """
    source = "{}_{}".format(db_name, db_version)
    with HashedSVIO(db_path) as reader:
        for record in reader:
            try:
                up_gene = selectAnnotSymbol(record[up_title], annotation_symbols, aliases_by_symbol)
                down_gene = selectAnnotSymbol(record[down_title], annotation_symbols, aliases_by_symbol)
                fusion_partners = "{}_@_{}".format(up_gene, down_gene)
                if fusion_partners not in fusions_by_partners:
                    fusions_by_partners[fusion_partners] = set()
                fusions_by_partners[fusion_partners].add(source)
            except Exception:
                print("warning", db_name, db_version, record["up_gene"], record["down_gene"], sep="\t")


def loadMitelman(db_path, db_version, fusions_by_partners, aliases_by_symbol, annotation_symbols):
    """
    Set fusions partners data from Mitelman database: MBCA.TXT.DATA,REF.TXT.DATA.

    :param db_path: Path to the Mitelman database MBCA.TXT.DATA,REF.TXT.DATA (format: TSV).
    :type db_path: str
    :param db_version: Database version to traceback sources in fusions_by_partners.
    :type db_version: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (chimerakb, cosmic, mitelman, pubmed, ...).
    :type fusions_by_partners: dict
    :param aliases_by_symbol: Gene name aliases by symbol.
    :type aliases_by_symbol: dict
    :param annotation_symbols: List of genes names known in genes annotations file.
    :type annotation_symbols: set
    """
    mbca_path = db_path
    pubmed_by_fusion = {}
    if "," in db_path:
        mbca_path, ref_path = db_path.split(",")
        pubmed_by_fusion = pubmedByFusion(ref_path)
    # MolClin    RefNo    InvNo    Morph    Topo    Immunology    GeneLength    GeneShort    GeneLong    KaryLength    KaryShort    KaryLong
    with HashedSVIO(mbca_path) as reader:
        for record in reader:
            if record["GeneShort"] != "":
                for fusion in record["GeneShort"].split(","):
                    if "/" in fusion:
                        try:
                            genes = fusion.replace("+", "").split("/")  # PDRG1/ARF3/RUNX1 => fusion between 3 genes
                            for up_gene, down_gene in zip(genes, genes[1:]):  # For each breakpoint
                                up_gene = selectAnnotSymbol(up_gene, annotation_symbols, aliases_by_symbol)
                                down_gene = selectAnnotSymbol(down_gene, annotation_symbols, aliases_by_symbol)
                                fusion_partners = "{}_@_{}".format(up_gene, down_gene)
                                source = "mitelman_{}".format(db_version)
                                if fusion_partners not in fusions_by_partners:
                                    fusions_by_partners[fusion_partners] = {source: set()}
                                if source not in fusions_by_partners[fusion_partners]:
                                    fusions_by_partners[fusion_partners][source] = set()
                                fusions_by_partners[fusion_partners][source].add(int(record["RefNo"]))
                                if "PMID" not in fusions_by_partners[fusion_partners]:
                                    fusions_by_partners[fusion_partners]["PMID"] = set()
                                    for pmid in pubmed_by_fusion[record["RefNo"]]:
                                        fusions_by_partners[fusion_partners]["PMID"].add(int(pmid))
                        except Exception:
                            print("warning", "mitelman", db_version, "\t".join(fusion.split("/")), sep="\t")


def pubmedByFusion(in_ref):
    """
    Return Pubmed IDs by fusion partners from REF.TXT.DATA from Mitelman database.

    :param in_ref: Path to the REF.TXT.DATA from Mitelman database (format: TSV).
    :type in_ref: str
    :return: Pubmed IDs by fusion partners.
    :rtype: dict
    """
    # RefNo    TitleLength    TitleShort    TitleLong    Volume    Year    Journal    Text    Abbreviation    AuthorsLength    AuthorsShort    AuthorsLong    Flag    Pubmed
    pubmed_by_fusion = {}
    with HashedSVIO(in_ref) as reader:
        for record in reader:
            if record["RefNo"] not in pubmed_by_fusion:
                pubmed_by_fusion[record["RefNo"]] = set()
            if record["Pubmed"] != "":
                pubmed_by_fusion[record["RefNo"]].add(record["Pubmed"])
    return pubmed_by_fusion


def writePartnersDb(db_path, fusions_by_partners):
    """
    Write known fusions partners database.

    :param db_path: Path to the fusions partners database (format: TSV).
    :type db_path: str
    :param fusions_by_partners: By partners (upGene_@_downGene) the ids of fusions by source (chimerakb, cosmic, mitelman, pubmed, ...).
    :type fusions_by_partners: dict
    """
    with HashedSVIO(db_path, "w") as writer:
        writer.titles = ["5prim_gene", "3prim_gene", "sources"]
        for partners, entries_by_src in fusions_by_partners.items():
            up_gene, down_gene = partners.split("_@_")
            sources = [src + ":" + ",".join([str(elt) for elt in sorted(ids)]) for src, ids in entries_by_src.items()]
            writer.write({
                "5prim_gene": up_gene,
                "3prim_gene": down_gene,
                "sources": "|".join(sources)
            })


class InputsDatabases(argparse.Action):
    """Manages inputs-databases parameters."""

    def __call__(self, parser, namespace, values, option_string=None):
        fct_by_model = {
            "babiceanu": loadBabiceanu,
            "bodymap": loadBodyMap,
            "cosmic": loadCosmic,
            "chimerdb": loadChimerdb,
            "mitelman": loadMitelman
        }
        databases = []
        for db_arg in values:
            if db_arg.count(":") != 2:
                raise argparse.ArgumentTypeError(
                    'Argument "{}" is invalid. The format must be: "MODEL:VERSION:PATH".'.format(
                        db_arg
                    )
                )
            model, version, path = db_arg.split(":")
            if model.lower() not in fct_by_model:
                raise argparse.ArgumentTypeError(
                    'Database model "{}" is invalid. It must be selected in {}.'.format(
                        model, sorted(fct_by_model.keys())
                    )
                )
            databases.append({
                "model": model,
                "parser": fct_by_model[model.lower()],
                "path": path,
                "version": version
            })
        setattr(namespace, self.dest, databases)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Write unique known fusions partners database from multiple databases.')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-l', '--input-aliases', required=True, help="Path to the file containing aliases between genes symbols (format: TSV). Each line contains a symbol (Gene name) with one of these aliases (Gene Synonym). This file can be obtain from Ensembl's biomart.")
    group_input.add_argument('-a', '--input-annotations', required=True, help='Path to the file containing the genes annotations used in analysis (format: GTF).')
    group_input.add_argument('-d', '--inputs-databases', action=InputsDatabases, nargs='+', required=True, help='Paths to databases format: "model:version:path". Model must be in ["babiceanu", "bodymap", "chimerdb", "cosmic", "mitelman"]. Example: -d cosmic:91:~/cosmic91_fusion.tsv chimerdb:Kb3.0:~/chimerKb.tsv mitelman:2019:~/MBCA.TXT.DATA,~/REF.TXT.DATA.')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-o', '--output-database', required=True, help='Path to the fusions partners database (format: TSV).')
    args = parser.parse_args()

    # Logger
    logging.basicConfig(format='%(asctime)s -- [%(filename)s][pid:%(process)d][%(levelname)s] -- %(message)s')
    log = logging.getLogger(os.path.basename(__file__))
    log.setLevel(logging.INFO)
    log.info("Command: " + " ".join(sys.argv))

    # Process
    log.info("Load required symbols and aliases")
    annotation_symbols = annotSymbols(args.input_annotations)
    aliases_by_symbol = aliasesBySymbols(args.input_aliases)
    fusions_by_partners = {}
    for db in args.inputs_databases:
        log.info("Load {}:{}".format(db["model"], db["version"]))
        db["parser"](db["path"], db["version"], fusions_by_partners, aliases_by_symbol, annotation_symbols)
    writePartnersDb(args.output_database, fusions_by_partners)
    log.info("End of job")
