#!/usr/bin/env python3

__author__ = 'Frederic Escudie'
__copyright__ = 'Copyright (C) 2017 CHU Toulouse'
__license__ = 'GNU General Public License'
__version__ = '1.4.2'

import os
import sys
import time
import argparse
import subprocess
from subprocess import Popen, PIPE
import uuid

BIN_DIR = os.path.dirname(__file__)
os.environ['PATH'] = os.environ['PATH'] + os.pathsep + BIN_DIR


########################################################################
#
# FUNCTIONS
#
########################################################################
class Cmd:
    """
    Command wrapper.

    :copyright: FROGS's team INRA.
    """

    def __init__(self, program, description, exec_parameters, version_parameters=None, interpreter=None):
        """
        :param exec_parameters: The parameters to execute the program. Two possibles syntaxes.
        If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
        Otherwise the parameters will be added after the program in command line.
        :type exec_parameters: str
        :param version_parameters: The parameters to get the program version. Two possibles syntaxes.
        If the parameter contains the string '##PROGRAM##', this tag will be replaced by the program parameter before submit.
        Otherwise the parameters will be added after the program in command line.
        :type version_parameters: str
        :param interpreter: The specific interpreter used to call program. For example '/home/user/venv/bin/python3' for a python script with
        dependencies installed in the virtual environment /home/user/venv.
        :type interpreter: str
        """
        self.program = program
        self.description = description
        self.exec_parameters = exec_parameters
        self.version_parameters = version_parameters
        self.interpreter = interpreter

    def get_cmd(self):
        """
        Return the command line.

        :return: The command line.
        :rtype: str
        """
        exec_call = self.program
        if self.interpreter is not None:
            exec_call = self.interpreter + " " + exec_call
        cmd = None
        if '##PROGRAM##' in self.exec_parameters:
            cmd = self.exec_parameters.replace('##PROGRAM##', exec_call)
        else:
            cmd = exec_call + ' ' + self.exec_parameters
        return cmd

    def get_version(self, location='stdout'):
        """
        Return the program version number.

        :param location: If the version command returns the version number on 'stdout' or on 'stderr'.
        :type location: str
        :return: Version number if this is possible, otherwise this method return 'unknown'.
        :rtype: str
        """
        if self.version_parameters is None:
            return "unknown"
        else:
            try:
                cmd = self.program + ' ' + self.version_parameters
                if self.interpreter is not None:
                    cmd = self.interpreter + ' ' + cmd
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                if location == 'stderr':
                    return stderr.decode('ascii').strip()
                else:
                    return stdout.decode('ascii').strip()
            except Exception:
                raise Exception("Version cannot be retrieve for the software '" + self.program + "'.")

    def parser(self, log_file):
        """
        Parse the command results to add information in log_file.

        :param log_file: Path to the sample process log file.
        :type log_file: str
        """
        pass

    def submit(self, log_file=None):
        """
        Launch command, trace this action in log and parse results.

        :param log_file: Path to the sample process log file.
        :type log_file: str
        """
        # Log
        if log_file is not None:
            FH_log = Logger(log_file)
            FH_log.write('# ' + self.description + '\n')
            FH_log.write('\tSoftware:\n\t\t' + os.path.basename(self.program) + ' version: ' + self.get_version() + '\n')
            FH_log.write('\tCommand:\n\t\t' + self.get_cmd() + '\n')
            FH_log.write('\tExecution:\n\t\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n')
            FH_log.close()
        # Process
        subprocess.check_output(self.get_cmd(), shell=True)
        # Log
        if log_file is not None:
            FH_log = Logger(log_file)
            FH_log.write('\t\tend:   ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n')
            FH_log.close()
            # Post-process results
            self.parser(log_file)


class Logger:
    """
    Log file handler.

    :copyright: FROGS's team INRA.
    """

    def __init__(self, filepath=None):
        """
        :param filepath: The log filepath. [default : STDOUT]
        :type filepath: str
        """
        self.filepath = filepath
        self.file_handle = None
        if self.filepath is not None and self.filepath is not sys.stdout:
            self.file_handle = open(self.filepath, "a")
        else:
            self.file_handle = sys.stdout

    def __del__(self):
        """Close file handler when the logger is detroyed."""
        self.close()

    def close(self):
        """Close file handler."""
        if self.filepath is not None and self.filepath is not sys.stdout:
            if self.file_handle is not None:
                self.file_handle.close()
                self.file_handle = None

    def write(self, msg):
        """
        Write msg on file.

        :param msg: The message to write.
        :type msg: str
        """
        self.file_handle.write(msg)

    @staticmethod
    def static_write(filepath, msg):
        """
        Write msg on file.

        :param filepath: The log filepath. [default : STDOUT]
        :type filepath: str
        :param msg: The message to write.
        :type msg: str
        """
        if filepath is not None and filepath is not sys.stdout:
            FH_log = open(filepath, "a")
            FH_log.write(msg)
            FH_log.close()
        else:
            sys.stdout.write(msg)


class TmpFiles:
    """
    Manage temporary files.

    :copyright: FROGS's team INRA.
    :note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add("toto.fasta")
            ...
            tmp_log = tmpFiles.add("log.txt")
            ...
        finaly:
            tmpFiles.deleteAll()
    """

    def __init__(self, tmp_dir, prefix=None):
        """
        :param tmp_dir: The temporary directory path.
        :type tmp_dir: str
        :param prefix: The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        :type prefix: str
        """
        if prefix is None:
            prefix = str(uuid.uuid4())
        self.files = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        Add a temporary file.

        :param filename: The filename without prefix.
        :type filename: str
        :param prefix: The prefix added [default: TmpFiles.prefix].
        :type prefix: str
        :param dir: The directory path [default: TmpFiles.tmp_dir].
        :type dir: str
        :return: The filepath.
        :rtype: str
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        Delete the specified temporary file.

        :param filepath: The file path to delete.
        :type filepath: str
        """
        self.files.remove(filepath)
        if os.path.exists(filepath):
            os.remove(filepath)

    def deleteAll(self):
        """Delete all temporary files."""
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)


class SamtoolsIndex(Cmd):
    """Index alignment."""

    def __init__(self, in_aln):
        """
        :param in_aln: Path to the alignment file (format: BAM).
        :type in_aln: str
        """
        cmd_param = "" + \
            " index" + \
            " " + in_aln

        Cmd.__init__(self,
                     "samtools",
                     "Index alignment.",
                     cmd_param,
                     "--version")

    def get_version(self):
        """
        Return the program version number.

        :return: version number if this is possible, otherwise this method return 'unknown'.
        :rtype: str
        """
        output = Cmd.get_version(self)
        return output.split("\n")[0].split(" ")[1].strip()


class SplitBAMByRG(Cmd):
    """Split BAM by groups of non-overlapping amplicons."""

    def __init__(self, in_design, in_aln, out_pattern):
        """
        :param in_design: Path to the amplicons description file (format: BED).
        :type in_design: str
        :param in_aln: Path to the alignment file (format: BAM).
        :type in_aln: str
        :param out_pattern: The path pattern for the outputted alignments files (format: BAM). In this path the keyword "{GP}" is replace by the group name for each group.
        :type out_pattern: str
        """
        cmd_param = "" + \
            " --input-design " + in_design + \
            " --input-aln " + in_aln + \
            " --output-pattern " + out_pattern

        Cmd.__init__(self,
                     "splitBAMByRG.py",
                     "Splits BAM by groups of non-overlapping amplicons.",
                     cmd_param,
                     "--version")


class AddRGOnBAM(Cmd):
    """Add tag on reads by origin (amplicon ID)."""

    def __init__(self, in_aln, out_aln, platform, sample, library):
        """
        :param in_aln: Path to the alignments file (format: BAM).
        :type in_aln: str
        :param out_aln: Path to the outputted alignments file (format: BAM).
        :type out_aln: str
        :param platform: Platform/technology used to produce the reads.
        :type platform: str
        :param sample: Sample. Use pool name where a pool is being sequenced.
        :type sample: str
        :param library: Library.
        :type library: str
        """
        cmd_param = "" + \
            " --pl " + platform + \
            " --sm " + sample + \
            " --lb " + library + \
            " --input-aln " + in_aln + \
            " --output-aln " + out_aln

        Cmd.__init__(self,
                     "addRGOnBAM.py",
                     "Adds tag on reads by origin (amplicon ID).",
                     cmd_param,
                     "--version")


class VarDictStep1(Cmd):
    """Dicover variants."""

    def __init__(self, in_reference, in_regions, in_aln, out_file, min_alt_freq=0.02, min_alt_count=4, min_base_qual=25, vardict_call="vardict-java"):
        """
        :param in_reference: Path to the reference sequences file (format: fasta).
        :type in_reference: str
        :param in_regions: Path to the amplicons design (format: BED). Start and end of the amplicons must be with primers.
        :type in_regions: str
        :param in_aln: Path to the alignments file (format: BAM).
        :type in_aln: str
        :param out_file: Path to the outputted file.
        :type out_file: str
        :param min_alt_freq: The threshold for allele frequency.
        :type min_alt_freq: float
        :param min_alt_count: The threshold for allele count.
        :type min_alt_count: int
        :param min_base_qual: The phred score for a base to be considered a good call.
        :type min_base_qual: int
        :param vardict_call: Command used to call vardict.
        :type vardict_call: str
        """
        cmd_param = "" + \
            " -r " + str(min_alt_count) + \
            " -f " + str(min_alt_freq) + \
            " -q " + str(min_base_qual) + \
            " -P 0" + \
            " -F 0" + \
            " -c 1 -S 2 -E 3 -g 4" + \
            " -b " + in_aln + \
            " -G " + in_reference + \
            " " + in_regions + \
            " > " + out_file

        Cmd.__init__(self,
                     vardict_call,
                     "Dicovers variants.",
                     cmd_param)


class VarDictStep2(Cmd):
    """Filter variant on strand bias."""

    def __init__(self, in_file, out_file):
        """
        :param in_file: Path to the input file.
        :type in_file: str
        :param out_file: Path to the outputted file.
        :type out_file: str
        """
        cmd_param = "" + \
            " cat " + in_file + " | " + \
            " ##PROGRAM##" + \
            " > " + out_file

        Cmd.__init__(self,
                     "teststrandbias.R",
                     "Filters variants on strand bias.",
                     cmd_param)


class VarDictStep3(Cmd):
    """Filter variants and converts to VCF."""

    def __init__(self, in_file, out_variants, min_alt_freq=0.02):
        """
        :param in_file: Path to the input file.
        :type in_file: str
        :param out_variants: Path to the outputted file (format: VCF).
        :type out_variants: str
        :param min_alt_freq: The threshold for allele frequency.
        :type min_alt_freq: float
        """
        cmd_param = "" + \
            " cat " + in_file + " | " + \
            " ##PROGRAM##" + \
            " -A" + \
            " -a" + \
            " -E" + \
            " -P 0" + \
            " -f " + str(min_alt_freq) + \
            " > " + out_variants

        Cmd.__init__(self,
                     "var2vcf_valid.pl",
                     "Filters variants and converts to VCF.",
                     cmd_param)


class GatherOverlappingRegions(Cmd):
    """Gather variants from non-overlapping groups."""

    def __init__(self, in_regions, in_variants, in_aln, out_variants):
        """
        :param in_regions: Path to the amplicons design. Start and end of the amplicons must be without primers (format: BED).
        :type in_regions: str
        :param in_variants: Path to the variants files (format: VCF).
        :type in_variants: str
        :param in_aln: Path to the alignments files (format: BAM). Each alignment file correspond to a variants file.
        :type in_aln: str
        :param out_variants: Path to the outputted file (format: VCF).
        :type out_variants: str
        """
        cmd_param = "" + \
            " --input-designs " + " ".join(in_regions) + \
            " --input-variants " + " ".join(in_variants) + \
            " --input-aln " + " ".join(in_aln) + \
            " --output-variants " + out_variants

        Cmd.__init__(self,
                     "mergeVCFAmpli.py",
                     "Gathers variants from non-overlapping groups.",
                     cmd_param,
                     "--version")


class MeltOverlappingRegions(Cmd):
    """Melt all the samples contained in variant file in one sample."""

    def __init__(self, spl_name, in_variants, out_variants):
        """
        :param spl_name: Name of the final sample.
        :type spl_name: str
        :param in_variants: Path to the variants file with several samples (format: VCF).
        :type in_variants: str
        :param out_variants: Path to the outputted file (format: VCF).
        :type out_variants: str
        """
        cmd_param = "" + \
            " --new-spl-name '" + spl_name + "'" + \
            " --input-variants " + in_variants + \
            " --output-variants " + out_variants

        Cmd.__init__(self,
                     "meltVCFSamples.py",
                     "Melts all the samples contained in variant file in one sample.",
                     cmd_param,
                     "--version")


class FilterVCFPrimers(Cmd):
    """Remove variants located on amplicons primers."""

    def __init__(self, in_sequences, in_regions, in_variants, out_variants):
        """
        :param in_sequences: Path to the reference sequences file (format: fasta). The reference used to discover variants.
        :type in_sequences: str
        :param in_regions: Path to the amplicons design with their primers (format: BED). The zone of interest is defined by thickStart and thickEnd. The amplicons must not have any overlap between them.
        :type in_regions: str
        :param in_variants: Path to the variants file (format: VCF). This file should be sorted by coordinates otherwise the execution time will be dramatically increased.
        :type in_variants: str
        :param out_variants: Path to the outputted variants file (format: VCF).
        :type out_variants: str
        """
        cmd_param = "" + \
            " --input-sequences " + in_sequences + \
            " --input-regions " + in_regions + \
            " --input-variants " + in_variants + \
            " --output-variants " + out_variants

        Cmd.__init__(self,
                     "filterVCFPrimers.py",
                     "Removes variants located on amplicons primers.",
                     cmd_param,
                     "--version")


class fixVardictHeader(Cmd):
    """Fix errors in VCF format in VarDict outputs."""

    def __init__(self, in_variants, out_variants):
        """
        :param in_variants: Path to the variants file (format: VCF).
        :type in_variants: str
        :param out_variants: Path to the outputted variants file (format: VCF).
        :type out_variants: str
        """
        cmd_param = "" + \
            " --variant-caller vardict" + \
            " --input-variants " + in_variants + \
            " --output-variants " + out_variants

        Cmd.__init__(self,
                     "fixVCallerVCF.py",
                     "Fix errors in VCF format in VarDict outputs.",
                     cmd_param,
                     "--version")


def filterBED(in_bed, in_names, out_bed, nb_col=None):
    """
    Filter a BED file with the list of names of regions to keep.

    :param in_bed: Path to the initial file (format: BED).
    :type in_bed: str
    :param in_names: Path to the file containing the list of names of the retained regions.
    :type in_names: str
    :param out_bed: Path to the filtered file (format: BED).
    :type out_bed: str
    :param nb_col: Number of columns in output.
    :type nb_col: int
    """
    # Retrieve retained regions names
    retained_regions = dict()
    with open(in_names) as FH_names:
        for line in FH_names:
            retained_regions[line.strip()] = 1
    # Filter BED
    with open(in_bed) as FH_in:
        with open(out_bed, "w") as FH_out:
            for line in FH_in:
                if line.startswith("browser ") or line.startswith("track ") or line.startswith("#"):
                    FH_out.write(line)
                else:
                    fields = [field.strip() for field in line.split("\t")]
                    if fields[3] in retained_regions:
                        if nb_col is not None:
                            line = "\t".join(fields[:nb_col]) + "\n"
                        FH_out.write(line)


def VarDictFct(in_reference, in_regions, in_aln, out_variants, logger, tmp_file, min_alt_freq=0.02, min_alt_count=4, min_base_qual=25, vardict_call="vardict-java"):
    """
    Dicover amplicons variants with VarDict.

    :param in_reference: Path to the reference sequences file (format: fasta).
    :type in_reference: str
    :param in_regions: Path to the amplicons design (format: BED). Start and end of the amplicons must be with primers.
    :type in_regions: str
    :param in_aln: Path to the alignments file (format: BAM).
    :type in_aln: str
    :param out_variants: Path to the outputted file (format: VCF).
    :type out_variants: str
    :param logger: Logger used to trace sub-commands.
    :type logger: Logger
    :param tmp_file: Temporaries files manager.
    :type tmp_file: TmpFiles
    :param min_alt_freq: The threshold for allele frequency.
    :type min_alt_freq: float
    :param min_base_qual: The phred score for a base to be considered a good call.
    :type min_base_qual: int
    :param vardict_call: Command used to call vardict.
    :type vardict_call: str
    """
    out_vardict = tmp.add("vardict.txt")
    out_strand_bias = tmp.add("strdBias.txt")
    VarDictStep1(in_reference, in_regions, in_aln, out_vardict, min_alt_freq, min_alt_count, min_base_qual, vardict_call).submit(logger)
    VarDictStep2(out_vardict, out_strand_bias).submit(logger)
    VarDictStep3(out_strand_bias, out_variants, min_alt_freq).submit(logger)


########################################################################
#
# MAIN
#
########################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description='Varaint calling on Illumina amplicon sequencing. It use VarDictJava (see: https://github.com/AstraZeneca-NGS/VarDictJava).')
    parser.add_argument('-m', '--min-alt-freq', default=0.02, type=float, help='Variants with an allele frequency under this value are not emitted. [Default: %(default)s]')
    parser.add_argument('-c', '--min-alt-count', default=4, type=int, help='Variants with an allele count under this value are not emitted. [Default: %(default)s]')
    parser.add_argument('-q', '--min-base-qual', default=25, type=int, help='The phred score for a base to be considered a good call. [Default: %(default)s]')
    parser.add_argument('-t', '--vardict-call', default="vardict-java", help='Command used to call vardict. [Default: %(default)s]')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    group_reference = parser.add_argument_group('Reference')  # Reference
    group_reference.add_argument('-g', '--input-genome', required=True, help='The path to the reference genome (format: fasta).')
    group_panel = parser.add_argument_group('Design')  # Design
    group_panel.add_argument('-pi', '--input-design-with-primers', required=True, help='The path to the amplicons design with their primers (format: BED).')
    group_panel.add_argument('-po', '--input-design-wout-primers', required=True, help='The path to the amplicons design without their primers (format: BED).')
    group_panel.add_argument('-pg', '--input-non-overlapping-design', required=True, help='The path to the list of amplicons (format: TSV). The first column is the ID of the amplicon and the second is the name of the group where the amplicon has no overlap with other amplicons of this group.')
    group_input = parser.add_argument_group('Inputs')  # Inputs
    group_input.add_argument('-l', '--library-name', help='The library name (example: "patient1_libA").')
    group_input.add_argument('-a', '--input-aln', required=True, help='The path to the alignment file (format: BAM).')
    group_output = parser.add_argument_group('Outputs')  # Outputs
    group_output.add_argument('-ov', '--output-variants', required=True, help='The path to the outputted variants file (format: VCF).')
    group_output.add_argument('-ol', '--output-log', default=sys.stdout, help='The path to the outputted log file (format: txt). [Default: STDOUT]')
    args = parser.parse_args()

    Logger.static_write(args.output_log, "## Application\n\tSoftware:\n\t\t" + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\n\tCommand:\n\t\t" + " ".join(sys.argv) + "\n\n")
    tmp = TmpFiles(os.path.dirname(args.output_variants))
    library_name = os.path.basename(args.input_aln).split(".")[0] if args.library_name is None else args.library_name

    # Get non-overlapping groups
    groups_names = set()
    with open(args.input_non_overlapping_design) as FH_gp:
        for line in FH_gp:
            if not line.startswith("#"):
                amplicon_id, group_name = [elt.strip() for elt in line.split("\t")]
                groups_names.add(group_name)
    groups_names = list(groups_names)

    # Split BAM in non-overlapping regions
    gp_alignment = [tmp.add(gp + ".bam") for gp in groups_names]
    out_bam_pattern = gp_alignment[-1][:-(len(groups_names[-1]) + 4)] + "{GP}.bam"
    SplitBAMByRG(args.input_non_overlapping_design, args.input_aln, out_bam_pattern).submit(args.output_log)

    # Variant calling
    groups = list()
    for idx_gp, curr_gp in enumerate(groups_names):
        curr_gp_aln = gp_alignment[idx_gp]

        # Index BAM
        tmp.files.append(curr_gp_aln + ".bai")
        SamtoolsIndex(curr_gp_aln).submit(args.output_log)

        # Select regions of current group
        curr_gp_regions = tmp.add(curr_gp + "_amplicons.txt")
        subprocess.check_output('grep "' + curr_gp + '$" ' + args.input_non_overlapping_design + ' | cut -f 1 > ' + curr_gp_regions, shell=True)
        curr_gp_regions_with_prim = tmp.add(curr_gp + "_withPrimers.bed")
        filterBED(args.input_design_with_primers, curr_gp_regions, curr_gp_regions_with_prim)
        curr_gp_regions_with_prim_4_col = tmp.add(curr_gp + "_withPrimers_4col.bed")
        filterBED(args.input_design_with_primers, curr_gp_regions, curr_gp_regions_with_prim_4_col, 4)
        curr_gp_regions_wout_prim = tmp.add(curr_gp + "_woutPrimers.bed")
        filterBED(args.input_design_wout_primers, curr_gp_regions, curr_gp_regions_wout_prim)

        # Add RG on BAM
        curr_gp_aln_new_RG = tmp.add(curr_gp + "_RG.bam")
        AddRGOnBAM(curr_gp_aln, curr_gp_aln_new_RG, "ILLUMINA", library_name, library_name).submit(args.output_log)
        tmp.files.append(curr_gp_aln_new_RG + ".bai")
        SamtoolsIndex(curr_gp_aln_new_RG).submit(args.output_log)

        # Call variants
        curr_gp_vcf = tmp.add(curr_gp + ".vcf")
        VarDictFct(
            args.input_genome,
            curr_gp_regions_with_prim_4_col,
            curr_gp_aln_new_RG, curr_gp_vcf,
            args.output_log,
            tmp,
            args.min_alt_freq,
            args.min_alt_count,
            args.min_base_qual,
            args.vardict_call
        )

        # Filters variants located on primers
        curr_gp_clean_vcf = tmp.add(curr_gp + "_clean.vcf")
        FilterVCFPrimers(args.input_genome, curr_gp_regions_with_prim, curr_gp_vcf, curr_gp_clean_vcf).submit(args.output_log)

        # Store current group files
        groups.append({
            "name": curr_gp,
            "aln": curr_gp_aln,
            "vcf": curr_gp_clean_vcf,
            "design_wout_prim": curr_gp_regions_wout_prim,
            "design_with_prim": curr_gp_regions_with_prim
        })

    # Merge overlapping amplicons
    out_gather = tmp.add("gatherOverlapping.vcf")
    GatherOverlappingRegions(
        [curr_gp["design_wout_prim"] for curr_gp in groups],
        [curr_gp["vcf"] for curr_gp in groups],
        [curr_gp["aln"] for curr_gp in groups],
        out_gather
    ).submit(args.output_log)
    out_melt = tmp.add("melt.vcf")
    MeltOverlappingRegions(library_name, out_gather, out_melt).submit(args.output_log)

    # Fix vardict header
    out_gather = tmp.add("gatherOverlapping.vcf")
    fixVardictHeader(out_melt, args.output_variants).submit(args.output_log)

    # Clean temporary files
    tmp.deleteAll()
