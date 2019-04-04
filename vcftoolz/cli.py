#!/usr/bin/env python

"""Console script for VCF Toolz."""

from __future__ import print_function
from __future__ import absolute_import

import argparse
import logging
import sys

from vcftoolz import vcftoolz
from vcftoolz.__init__ import __version__

# Ignore flake8 errors in this module
# fl ake8: noqa


def parse_arguments(system_args):
    """Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    description = "Tools for working with Variant Call Format files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--version", action="version", version="%(prog)s version " + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    help_str = "Compare VCF files."
    description = """Compare and analyze the variants found in multiple input VCF files. A summary report of differences is printed to
      stdout.  Venn diagrams in pdf format are generated in the current working directory for up to 6 VCF files. A spreadsheet of
      tabulated results can be optionally generated.  When the truth is known, the summary file contains additional metrics describing
      the effectiveness of the variant caller as well as the effectiveness of the filters."""
    subparser = subparsers.add_parser("compare", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(nargs=1,              dest="vcf_path1",        type=str, metavar="VCF",  help="VCF file")
    subparser.add_argument(nargs='+',            dest="vcf_path_list",    type=str, metavar="VCF",  help="VCF file")
    subparser.add_argument("--exclude_snps",     dest="exclude_snps",     action="store_true",      help="Exclude snp calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_indels",   dest="exclude_indels",   action="store_true",      help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_vars",     dest="exclude_vars",     action="store_true",      help="Exclude variants other than snps and indels.")
    subparser.add_argument("--exclude_refs",     dest="exclude_refs",     action="store_true",      help="Exclude reference calls.")
    subparser.add_argument("--exclude_hetero",   dest="exclude_hetero",   action="store_true",      help="Exclude heterozygous calls.")
    subparser.add_argument("--exclude_filtered", dest="exclude_filtered", action="store_true",      help="Exclude filtered calls (FT or FILTER is not PASS).")
    subparser.add_argument("--exclude_missing",  dest="exclude_missing",  action="store_true",      help="Exclude calls with all data elements missing.")
    subparser.add_argument("--truth",            dest="truth_flag",       action="store_true",      help="Additional metrics are generated assuming the first VCF file is the truth. This also triggers extra analysis of filtered calls.")
    subparser.add_argument("--tabulate",         dest="table_file",       type=str, metavar='FILE', help="Tabulate the results in the specified tab-separated-value file.")
    subparser.set_defaults(func=compare_command)

    help_str = "Reformat VCF data into a tall, narrow format."
    description = "Convert a VCF file into a tall and narrow format.  Output is printed to stdout."
    subparser = subparsers.add_parser("narrow", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(                      dest="vcf_path",         type=str, metavar="VCF",  help="VCF file")  # noqa: E201
    subparser.add_argument("--exclude_snps",     dest="exclude_snps",     action="store_true",      help="Exclude snp calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_indels",   dest="exclude_indels",   action="store_true",      help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_vars",     dest="exclude_vars",     action="store_true",      help="Exclude variants other than snps and indels.")
    subparser.add_argument("--exclude_refs",     dest="exclude_refs",     action="store_true",      help="Exclude reference calls.")
    subparser.add_argument("--exclude_hetero",   dest="exclude_hetero",   action="store_true",      help="Exclude heterozygous calls.")
    subparser.add_argument("--exclude_filtered", dest="exclude_filtered", action="store_true",      help="Exclude filtered calls (FT or FILTER is not PASS).")
    subparser.add_argument("--exclude_missing",  dest="exclude_missing",  action="store_true",      help="Exclude calls with all data elements missing.")
    subparser.set_defaults(func=narrow_command)

    help_str = "Count the number of positions, calls, variants, filters, etc."
    description = "Count samples, positions, calls, snps, indels, other variants, filtered calls, missing calls, and filter reasons.  Output is printed to stdout."
    subparser = subparsers.add_parser("count", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(                      dest="vcf_path",         type=str, metavar="VCF",  help="VCF file")  # noqa: E201
    subparser.add_argument("--exclude_snps",     dest="exclude_snps",     action="store_true",      help="Exclude snp calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_indels",   dest="exclude_indels",   action="store_true",      help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_vars",     dest="exclude_vars",     action="store_true",      help="Exclude variants other than snps and indels.")
    subparser.add_argument("--exclude_refs",     dest="exclude_refs",     action="store_true",      help="Exclude reference calls.")
    subparser.add_argument("--exclude_hetero",   dest="exclude_hetero",   action="store_true",      help="Exclude heterozygous calls.")
    subparser.add_argument("--exclude_filtered", dest="exclude_filtered", action="store_true",      help="Exclude filtered calls (FT or FILTER is not PASS).")
    subparser.add_argument("--exclude_missing",  dest="exclude_missing",  action="store_true",      help="Exclude calls with all data elements missing.")
    subparser.set_defaults(func=count_command)

    help_str = "Plot calls along the length of the genome."
    description = "Plot calls along the length of the genome.  By default, all contigs are joined together in order of decreasing size.  This is intended to be used for pooled vcf files with many samples."
    subparser = subparsers.add_parser("plot", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(                      dest="vcf_path",         type=str, metavar="VCF",    help="VCF file")  # noqa: E201
    subparser.add_argument(                      dest="reference_path",   type=str, metavar="FASTA",  help="Reference fasta file")  # noqa: E201
    subparser.add_argument(                      dest="output_path",      type=str, metavar="OUTPUT", help="Output file.  You should use the extension .pdf or .png")  # noqa: E201
    subparser.add_argument("--chrom",            dest="chrom",            default=None,               help="Only plot a single named contig. By default all contigs are joined together in order of decreasing size.")
    subparser.add_argument("--exclude_snps",     dest="exclude_snps",     action="store_true",        help="Exclude snp calls. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_indels",   dest="exclude_indels",   action="store_true",        help="Exclude insertions and deletions. A heterozygous call with both snp and indel is not excluded unless both snps and indels are excluded.")
    subparser.add_argument("--exclude_vars",     dest="exclude_vars",     action="store_true",        help="Exclude variants other than snps and indels.")
    subparser.add_argument("--exclude_refs",     dest="exclude_refs",     action="store_true",        help="Exclude reference calls.")
    subparser.add_argument("--exclude_hetero",   dest="exclude_hetero",   action="store_true",        help="Exclude heterozygous calls.")
    subparser.add_argument("--exclude_filtered", dest="exclude_filtered", action="store_true",        help="Exclude filtered calls (FT or FILTER is not PASS).")
    subparser.add_argument("--exclude_missing",  dest="exclude_missing",  action="store_true",        help="Exclude calls with all data elements missing.")
    subparser.set_defaults(func=plot_command)

    args = parser.parse_args(system_args)
    return args


def compare_command(args):
    """Compare and analyze the snps found in two or more input VCF files.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    vcf_path_list = args.vcf_path1 + args.vcf_path_list
    if args.truth_flag and args.exclude_filtered:
        logging.info("The --exclude_filtered flag is ignored when --truth is set.  It will be handled automatically.  The --truth flag triggers extra analysis of filtered calls.")
    if args.truth_flag and len(vcf_path_list) != 2:
        logging.error("When the truth option is used, there must be exactly 2 VCF files.")
        exit(1)
    vcftoolz.compare(args.truth_flag, vcf_path_list, args.exclude_snps, args.exclude_indels, args.exclude_vars, args.exclude_refs, args.exclude_hetero, args.exclude_filtered, args.exclude_missing, args.table_file)


def narrow_command(args):
    """Convert a VCF file into a tab delimited set of snp calls, one per line.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    vcftoolz.narrow(args.vcf_path, args.exclude_snps, args.exclude_indels, args.exclude_vars, args.exclude_refs, args.exclude_hetero, args.exclude_filtered, args.exclude_missing)


def count_command(args):
    """Count the number of positions and calls.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    vcftoolz.count(args.vcf_path, args.exclude_snps, args.exclude_indels, args.exclude_vars, args.exclude_refs, args.exclude_hetero, args.exclude_filtered, args.exclude_missing)


def plot_command(args):
    """Plot calls along the length of the genome.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.
    """
    vcftoolz.plot(args.vcf_path, args.reference_path, args.output_path, args.exclude_snps, args.exclude_indels, args.exclude_vars, args.exclude_refs, args.exclude_hetero, args.exclude_filtered, args.exclude_missing, args.chrom)


def run_command_from_args(args):
    """Run a subcommand with previously parsed arguments in an argparse namespace.

    This function is intended to be used for unit testing.

    Parameters
    ----------
    args : Namespace
        Command line arguments are stored as attributes of a Namespace.
        The args are obtained by calling parse_argument_list().

    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def run_from_line(line):
    """Run a command with a command line.

    This function is intended to be used for unit testing.

    Parameters
    ----------
    line : str
        Command line.

    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    argv = line.split()
    args = parse_arguments(argv)
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def main():
    """This is the main function which is turned into an executable
    console script by the setuptools entry_points.  See setup.py.

    To run this function as a script, first install the package:
        $ python setup.py develop
        or
        $ pip install --user vcftoolz

    Parameters
    ----------
    This function must not take any parameters

    Returns
    -------
    The return value is passed to sys.exit().
    """
    enable_log_timestamps = False
    if enable_log_timestamps:
        logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO)
    else:
        logging.basicConfig(format="%(message)s", level=logging.INFO)
    args = parse_arguments(sys.argv[1:])
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


# This snippet lets you run the cli without installing the entrypoint.
if __name__ == "__main__":
    sys.exit(main())
