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
* Venn Diagram
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
import itertools
import collections

__version__ = '0.1.0'


def report_error(message):
    """
    Send an error message to stderr.
    """
    print(message, file=sys.stderr)


def get_unique_set_elements(sets):
    """
    Given a list of sets, return a new list of sets with the overlapping
    elements removed, leaving only the set elements that are unique to each
    set.

    Parameters
    ----------
    sets : list of set

    Returns
    -------
    unique_element_sets : list of set
    """
    # To avoid quadratic runtime, make an initial pass to determine which elements appear in only one set
    element_counts = collections.Counter(itertools.chain.from_iterable(sets))

    all_unique_elements = {element for element, count in element_counts.items() if count == 1}
    unique_element_sets = [element_set & all_unique_elements for element_set in sets]
    return unique_element_sets


def get_missing_set_elements(sets):
    """
    Given a list of sets, return a new list of sets containing the elements that are
    missing from each set but are present in more than one other set.

    Parameters
    ----------
    sets : list of set

    Returns
    -------
    missing_element_sets : list of set
    """
    # To avoid quadratic runtime, make an initial pass to determine which elements appear in more than one set
    element_counts = collections.Counter(itertools.chain.from_iterable(sets))

    all_duplicate_elements = {element for element, count in element_counts.items() if count > 1}
    missing_element_sets = [all_duplicate_elements - element_set for element_set in sets]
    return missing_element_sets



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

    Returns
    -------
    samples : list of str
        List of sample names in the VCF file
    """
    with open(vcf_path) as inp:
        reader = vcf.VCFReader(inp)
        return reader.samples


class SampleSnps(object):
    """
    Class to hold a list of snps and the genotype
    data per snp position.
    """

    def __init__(self, snp_list, alt_dict, format_dict, call_dict):
        """
        Initialize the Snps container.

        Parameters
        ----------
        snp_list : list of tuples
            List of tuples of (chrom, pos, ref, call, sample)
        alt_dict :  dict
            Dictionary mapping from keyed snp tuple to actual alt list
        format_dict : dict
            Dictionary mapping from keyed snp tuple to FORMAT string
        call_dict : dict
            Dictionary mapping from keyed snp tuple to namedtuple of genotype data elements
        """
        self.snp_list = snp_list
        self.alt_dict = alt_dict
        self.format_dict = format_dict
        self.call_dict = call_dict


def get_snp_list(vcf_path, all_records, passfilter):
    """
    Get the list of snps in a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    all_records : bool
        Flag to cause all processing of all records, regardless of whether the
        record is a snp record
    passfilter : bool
        Process only the VCF samples or records having PASS FT element or PASS
        filter.  The filter element is always ignored when samples have the FT
        element regardless of this option.

    Returns
    -------
    snps : SampleSnps
        Container of:
            snp_list : list of tuples of (chrom, pos, ref, alt, sample)
            alt_dict :  dictionary mapping from keyed snp tuple to actual alt list
            format_dict : dictionary mapping from keyed snp tuple to FORMAT string
            call_dict : dictionary mapping from keyed snp tuple to namedtuple of genotype data elements
    """
    with open(vcf_path) as inp:
        reader = vcf.VCFReader(inp)
        snps = []
        alt_dict = dict()
        format_dict = dict()
        call_dict = dict()
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

                # If there are filters per sample, use them
                try:
                    FT = sample.data.FT
                    if FT and FT != "PASS":
                        continue
                # Otherwise, use the filter for the whole record
                except:
                    if passfilter and len(record.FILTER) > 0:
                        continue

                snp = (record.CHROM, int(record.POS), record.REF, bases, sample.sample)
                snps.append(snp)
                alt_dict[snp] = record.ALT
                format_dict[snp] = record.FORMAT
                call_dict[snp] = sample.data
    return SampleSnps(snps, alt_dict, format_dict, call_dict)


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

    parser.add_argument(dest="vcf_path_list", type=str, metavar="VcfFile", nargs='+',  help="List of VCF files")
    parser.add_argument("-a", "--allrecords", action='store_true', dest="all_records", help="Process all VCF records assuming all records are snp records.")
    parser.add_argument("-p", "--pass",       action='store_true', dest="passfilter",  help="Process only the VCF samples or records having PASS FT element or PASS filter.  The filter element is always ignored when samples have the FT element regardless of this option.")
    parser.add_argument('--version',          action='version', version='%(prog)s version ' + __version__)

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
    num_vcf_files = len(base_vcf_file_name_list)

    # Get the list of sample names in each VCF file
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
    sample_snps_list = [get_snp_list(path, args.all_records, args.passfilter) for path in args.vcf_path_list] # List of SampleSnps
    snp_set_list = [set(sample_snps.snp_list) for sample_snps in sample_snps_list] # list of set of tuples of (chrom, pos, ref, call, sample)
    alt_dict_list = [sample_snps.alt_dict for sample_snps in sample_snps_list]
    format_dict_list = [sample_snps.format_dict for sample_snps in sample_snps_list]
    call_dict_list = [sample_snps.call_dict for sample_snps in sample_snps_list]

    # Find the snps that only appear in one of the VCF files
    unique_snps_sets = get_unique_set_elements(snp_set_list)

    # Find the snps that are missing from each VCF file but are present in more than one other VCF file
    if num_vcf_files >= 3:
        missing_snps_sets = get_missing_set_elements(snp_set_list)

    # Print some statistics
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(sample_set_list[i])
        print("{count}\tSamples in {dataset}".format(dataset=dataset, count=count))
    print()

    position_set_list = []
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        positions = set([t[1] for t in snp_set_list[i]])
        position_set_list.append(positions)
        count = len(positions)
        print("{count}\tPositions having snps in {dataset}".format(dataset=dataset, count=count))
    print()

    unique_positions_sets = get_unique_set_elements(position_set_list)
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(unique_positions_sets[i])
        print("{count}\tSNP Positions only in {dataset}".format(dataset=dataset, count=count))
    print()

    if num_vcf_files >= 3:
        missing_positions_sets = get_missing_set_elements(position_set_list)
        for i in range(num_vcf_files):
            dataset = base_vcf_file_name_list[i]
            count = len(missing_positions_sets[i])
            print("{count}\tSNP Positions missing in {dataset}, but present in at least 2 other VCF files".format(dataset=dataset, count=count))
        print()

    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(snp_set_list[i])
        print("{count}\tSample snps in {dataset}".format(dataset=dataset, count=count))
    print()

    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(unique_snps_sets[i])
        print("{count}\tSample snps only in {dataset}".format(dataset=dataset, count=count))
    print()

    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            dataset = base_vcf_file_name_list[i]
            count = len(missing_snps_sets[i])
            print("{count}\tSample snps missing in {dataset}, but present in at least 2 other VCF files".format(dataset=dataset, count=count))
        print()

    # Print the positions present in only one of the VCF files
    for i in range(num_vcf_files):
        print("\nPositions only in %s" % base_vcf_file_name_list[i])
        sorted_positions = sorted(list(unique_positions_sets[i]))
        if len(sorted_positions) == 0:
            print("None")
        else:
            for pos in sorted_positions:
                print(pos)

    # Print the positions missing in each of the VCF files
    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            print("\nPositions missing in %s, but present in at least 2 other VCF files:" % base_vcf_file_name_list[i])
            sorted_positions = sorted(list(missing_positions_sets[i]))
            if len(sorted_positions) == 0:
                print("None")
            else:
                for pos in sorted_positions:
                    print(pos)

    # Print the snps present in only one of the VCF files
    for i in range(num_vcf_files):
        print("\nSample snps only in %s" % base_vcf_file_name_list[i])
        sorted_snps = sorted(list(unique_snps_sets[i]))
        if len(sorted_snps) == 0:
            print("None")
        else:
            print("CHROM   \tPOS\tREF\tALT\tSAMPLE  \tFORMAT\tDATA")
            for snp in sorted_snps:
                fields = [str(x) for x in snp]
                fields[3] = ",".join(map(str, alt_dict_list[i][snp]))
                format_str = format_dict_list[i][snp]
                format_keys = format_str.split(":")
                call_data = call_dict_list[i][snp]
                call_data_str = ":".join([str(getattr(call_data, k, None)) for k in format_keys])
                fields.append(format_str)
                fields.append(call_data_str)
                print('\t'.join(fields))


    # Print the snps missing in each of the VCF files
    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            print("\nSample snps missing in %s, but present in at least 2 other VCF files:" % base_vcf_file_name_list[i])
            sorted_snps = sorted(list(missing_snps_sets[i]))
            if len(sorted_snps) == 0:
                print("None")
            else:
                print("CHROM   \tPOS\tREF\tALT\tSAMPLE")
                for snp in sorted_snps:
                    fields = [str(x) for x in snp]
                    print('\t'.join(fields))


    # Generate a venn diagram if the necessary packages are installed
    #generate_venn_diagrams(snp_set_list, base_vcf_file_name_list, 'snps.venn.pdf')
    if num_vcf_files == 2 or num_vcf_files == 3:
        try:
            from matplotlib import pyplot as plt
            from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
        except:
            report_error("Skipping venn diagram creation.  Requires matplotlib and matplotlib_venn.")
            exit(1)

        colors2 = ['red', 'blue', 'magenta']
        alpha2 = [0.6, 0.4, 0.1]
        if len(base_vcf_file_name_list) == 2:
            c = venn2(snp_set_list, set_labels=base_vcf_file_name_list)
            c.get_patch_by_id('10').set_color(colors2[0])
            c.get_patch_by_id('01').set_color(colors2[1])
            c.get_patch_by_id('11').set_color(colors2[2])
            c.get_patch_by_id('10').set_alpha(alpha2[0])
            c.get_patch_by_id('01').set_alpha(alpha2[1])
            c.get_patch_by_id('11').set_alpha(alpha2[2])
            plt.title("Venn Diagram of SNPs")
        else:
            draw_3_way = False
            fig, axes = plt.subplots(3 + draw_3_way)
            plt_idx = 0
            if draw_3_way:
                venn3(snp_set_list, set_labels=base_vcf_file_name_list, ax=axes[plt_idx])
                plt_idx += 1
            for pair in itertools.combinations(range(num_vcf_files), 2):
                sets = [snp_set_list[k] for k in range(num_vcf_files) if k in pair]
                names = [base_vcf_file_name_list[k] for k in range(num_vcf_files) if k in pair]
                c = venn2(sets, set_labels=names, ax=axes[plt_idx])
                c.get_patch_by_id('10').set_color(colors2[0])
                c.get_patch_by_id('01').set_color(colors2[1])
                c.get_patch_by_id('11').set_color(colors2[2])
                c.get_patch_by_id('10').set_alpha(alpha2[0])
                c.get_patch_by_id('01').set_alpha(alpha2[1])
                c.get_patch_by_id('11').set_alpha(alpha2[2])
                plt_idx += 1
            axes[0].set_title("Venn Diagram of SNPs")
        plt.show()
        plt.savefig('snps.venn.pdf')




if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    main(args)


