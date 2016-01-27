#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
Script to compare VCF files.

Inputs
------
two or more VCF files

Outputs
-------
* metrics per VCF file
* lists of snps exclusively found per VCF file
* Venn Diagram (planned for future release)
* parsimony trees, one per input VCF file (planned for future release)

Errors
------
* Sample names are not the same and are not compareable
* Different numnber of samples

Other Related Scripts
----------------------
compare-vcf.sh
    This script runs bcftools isec to compare two vcf files.
    it compresses and indexes the vcf files as needed for isec.

~/mnt/gnome3/DropBox/Steve/wgs-comparisions/analysis-tools/vcf-to-narrow.py
    Parses a VCF file to find snps

~/Projects/analyze-snpmutator-snppipeline-results/analysis.py
    Pandas code to compare and contrast data sets with inner and outer joins.

~/Projects/analyze-snpmutator-snppipeline-results/venn.100x.20x.fn.py
    Venn diagram example.

~/Projects/multi-vcf/plot-snps.2a.py
    Plot a histogram of SNPS along the length of 2 genomes assuming there are no overlapping snps.
"""

from __future__ import print_function
import argparse
import sys
import os
import vcf

__version__ = '0.1.0'


def report_error(message):
    """
    Send an error message to stderr.
    """
    print(message, file=sys.stderr)



def verify_non_empty_input_files(error_prefix, file_list):
    """
    Verify each file in a list of files exists and is non-empty.
    Missing or empty files are reported in the verbose log.

    Parameters
    ----------
    error_prefix : str
        first part of error message to be logged
    file_list : list 
        relative or absolute paths to files

    Returns
    -------
    bad_count : int
        number of missing or empty files
    """
    bad_count = 0
    for file_path in file_list:

        if not os.path.isfile(file_path):
            bad_count += 1
            err_message = "%s %s does not exist." % (error_prefix, file_path)
            report_error(err_message)
            continue
        if os.path.getsize(file_path) == 0:
            bad_count += 1
            err_message = "%s %s is empty." % (error_prefix, file_path)
            report_error(err_message)
            continue

    return bad_count


def get_sample_list(vcf_path):
    """
    Get the list of sample names in a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    """
    with open(vcf_path) as inp:
        reader = vcf.VCFReader(inp)
        return reader.samples


def get_snp_list(vcf_path, all_records):
    """
    Get the list of snps in a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    all_records : bool
        Flag to cause all processing of all records, regardless of whether the
        record is a snp record

    Returns
    -------
    snp_list : list of tuples
        List of tuples of (chrom, pos, ref, alt, sample)
    """
    with open(vcf_path) as inp:
        reader = vcf.VCFReader(inp)
        snps = []
        for record in reader:
            if not all_records and not record.is_snp:
                print("Warning: the vcf record is not a snp record at position %i in file %s" % (record.POS, vcf_path))
                #for a in record.ALT:
                #    print("ref=%s alt value=%s alt type=%s" % (record.REF, a, a.type))
                continue
            for sample in record.samples:
                if not sample.is_variant:
                    continue
                bases = sample.gt_bases
                if len(bases) == 3 and bases[0] == bases[2]: # e.g. G/G
                    bases = bases[0]
                if bases == "N":
                    continue
                try:
                    FT = sample.data.FT
                except:
                    FT = None
                if FT and FT != "PASS":
                    continue
                snp = (record.CHROM, int(record.POS), record.REF, bases, sample.sample)
                snps.append(snp)
    return snps


def parse_arguments(system_args):
    """
    Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    usage = """Compare and analyze the snps found in multiple input VCF files."""
               
    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(dest="vcf_path_list", type=str, metavar="VcfFile", nargs='+', help="List of VCF files")
    parser.add_argument("-a", "--allrecords", action='store_true', dest="all_records", help="Process all VCF records assuming all records are snp records.")
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)

    args = parser.parse_args(system_args)
    if len(args.vcf_path_list) < 2:
        parser.print_usage(sys.stderr)
        print("error: at least 2 files required", file=sys.stderr)
        sys.exit(1)
    return args


def main(args):
    """
    Compare and analyze the snps found in two or three input VCF files.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv, but can be set programmatically for unit testing
        or other purposes.

    See Also
    --------
    parse_arguments()
    """
    # Validate input files
    bad_files_count = verify_non_empty_input_files("VCF file", args.vcf_path_list)
    if bad_files_count > 0:
        exit(1)

    # Get short base file name for each input
    vcf_file_name_list = [os.path.basename(p) for p in args.vcf_path_list]
    base_vcf_file_name_list = [os.path.splitext(f)[0] for f in vcf_file_name_list]

    # Get the sample names
    sample_set_list = [set(get_sample_list(path)) for path in args.vcf_path_list]
    all_samples = set()
    all_samples = all_samples.union(*sample_set_list)
    missing_samples = [all_samples - s for s in sample_set_list]
    for i in range(len(base_vcf_file_name_list)):
        dataset = base_vcf_file_name_list[i]
        count = len(missing_samples[i])
        if count > 0:
            for sample in missing_samples[i]:
                print("Warning: dataset {dataset} is missing sample {sample}".format(dataset=dataset, sample=sample))
    if any(missing_samples):
        print()

    # Extract the snps from each vcf file
    snp_list_list = [get_snp_list(path, args.all_records) for path in args.vcf_path_list]
    snp_set_list = [set(snp_list) for snp_list in snp_list_list]
    unique_snps_sets = [u - s for u in snp_set_list for s in snp_set_list if not s is u]

    # Print some statistics
    for i in range(len(base_vcf_file_name_list)):
        dataset = base_vcf_file_name_list[i]
        count = len(sample_set_list[i])
        print("Number of samples in {dataset}:\t{count}".format(dataset=dataset, count=count))
    for i in range(len(base_vcf_file_name_list)):
        dataset = base_vcf_file_name_list[i]
        positions = set([t[1] for t in snp_set_list[i]])
        count = len(positions)
        print("Number of positions having snps in {dataset}:\t{count}".format(dataset=dataset, count=count))
    for i in range(len(base_vcf_file_name_list)):
        dataset = base_vcf_file_name_list[i]
        count = len(snp_set_list[i])
        print("Number of sample snps in {dataset}:\t{count}".format(dataset=dataset, count=count))
    for i in range(len(base_vcf_file_name_list)):
        dataset = base_vcf_file_name_list[i]
        count = len(unique_snps_sets[i])
        print("Number of sample snps only in {dataset}:\t{count}".format(dataset=dataset, count=count))


    # Print unique sets
    for i in range(len(base_vcf_file_name_list)):
        print("\nSample snps only in %s" % base_vcf_file_name_list[i])
        sorted_snps = sorted(list(unique_snps_sets[i]))
        if len(sorted_snps) == 0:
            print("None")
        else:
            print("CHROM   \tPOS\tREF\tALT\tSAMPLE")
            for snp in sorted_snps:
                stringfy_snp = [str(x) for x in snp]
                print('\t'.join(stringfy_snp))



if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)


