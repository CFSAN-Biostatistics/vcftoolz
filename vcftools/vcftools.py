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

Errors Reported
---------------
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

__version__ = '0.3.0'

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


def get_set_intersections(sets):
    """
    Given a list of sets, return a new list of sets with all the possible
    mutually exclusive overlapping combinations of those sets.  Another way
    to think of this is the mutually exclusive sections of a venn diagram
    of the sets.  If the original list has N sets, the returned list will
    have (2**N)-1 sets.

    Parameters
    ----------
    sets : list of set

    Returns
    -------
    combinations : list of tuple
        tag : str
            Binary string representing which sets are included / excluded in
            the combination.
        set : set
            The set formed by the overlapping input sets.
    """
    num_combinations = 2 ** len(sets)
    bit_flags = [2 ** n for n in range(len(sets))]
    flags_zip_sets = [z for z in zip(bit_flags, sets)]

    combo_sets = []
    for bits in range(num_combinations - 1, 0, -1):
        include_sets = [s for flag, s in flags_zip_sets if bits & flag]
        exclude_sets = [s for flag, s in flags_zip_sets if not bits & flag]
        combo = set.intersection(*include_sets)
        combo = set.difference(combo, *exclude_sets)
        tag = ''.join([str(int((bits & flag) > 0)) for flag in bit_flags])
        combo_sets.append((tag, combo))
    return combo_sets


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


SnpTuple = collections.namedtuple("SnpTuple", "chrom pos ref gt sample")


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
        snp_list : list of SnpTuple namedtuple
            List of SnpTuple namedtuple containing (chrom, pos, ref, gt, sample)
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
            snps : list of namedtuple SnpTuple containing (chrom, pos, ref, gt, sample)
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

                snp = SnpTuple(record.CHROM, int(record.POS), record.REF, bases, sample.sample)
                snps.append(snp)
                alt_dict[snp] = record.ALT
                format_dict[snp] = record.FORMAT
                call_dict[snp] = sample.data
    return SampleSnps(snps, alt_dict, format_dict, call_dict)


def tabulate_results(snp_sets, samples, table_file_path):
    """
    Create a TSV spreadsheet with one row per position and one column per sample
    showing the overlapping intersection of snps from each VCF file.

    Parameters
    ----------
    snp_sets : list
        list of set of SnpTuple namedtuples of (chrom, pos, ref, gt, sample).
        One set per VCF file in the same order as the VCF files on the command
        line.
    samples : list
        List of sample names
    table_file_path : str
        Path to the TSV file to be written.
    """
    with open(table_file_path, 'w') as f:
        # Remember the set of all the different genotypes found at sample positions where at least one VCF file contained a snp.
        # This is to handle the possibility where different genotypes are reported in different VCF files.
        pos_sample_gt = collections.defaultdict(set)
        for snps in snp_sets:
            for snp in snps:
                key = snp.chrom + str(snp.pos) + snp.sample
                pos_sample_gt[key].add(snp.gt)

        # Remember which VCF file contained each snp -- separate into tagged venn diagram sections
        pos_sample_venn_tag = dict()
        snp_set_list = [{(s.chrom, s.pos, s.sample) for s in snps} for snps in snp_sets]  # don't distinguish by gt
        snp_intersections = get_set_intersections(snp_set_list)
        for tag, snps in snp_intersections:
            for chrom, pos, sample in snps:
                key = chrom + str(pos) + sample
                pos_sample_venn_tag[key] = tag

        # Write the header row
        sorted_samples = sorted(samples)
        header_row = ["Chrom", "Pos", "Ref"] + sorted_samples
        f.write("%s\n" % '\t'.join(header_row))

        # Get the membership of each venn diagram
        position_set_list = [{(s.chrom, s.pos, s.ref) for s in snps} for snps in snp_sets]
        pos_intersections = get_set_intersections(position_set_list)
        for _, positions in pos_intersections:
            if len(positions) == 0:
                continue
            positions = sorted(positions) # set becomes sorted list
            for chrom, pos, ref in positions:
                # Prepare and write spreadsheet rows, one row per position and one column per sample
                pos_str = str(pos)
                spreadsheet_row = [chrom, pos_str, ref]

                # Lookup the intersection tag indicating which VCF files contained the snp
                spreadsheet_row_tag = [pos_sample_venn_tag.get(chrom + pos_str + sample, '') for sample in sorted_samples]

                # Lookup all the genotypes found at this sample position (usually just one or none at all)
                spreadsheet_row_gt = [','.join(sorted(pos_sample_gt.get(chrom + pos_str + sample, ['0']))) for sample in sorted_samples]

                # Join the tag to the gt list to make the spreadsheet row cells
                spreadsheet_row += [gt+'-'+tag if tag else gt for tag, gt in zip(spreadsheet_row_tag, spreadsheet_row_gt)]

                f.write("%s\n" % '\t'.join(spreadsheet_row))


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
    description = """VCF Tools."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version',          action='version', version='%(prog)s version ' + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand       ")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    description = "Compare and analyze the snps found in multiple input VCF files."
    subparser = subparsers.add_parser("compare", formatter_class=formatter_class, description=description, help="Compare VCF files.")
    subparser.add_argument(dest="vcf_path1",     type=str, metavar="VcfFile", nargs=1,  help="VCF file")
    subparser.add_argument(dest="vcf_path_list", type=str, metavar="VcfFile", nargs='+',  help="VCF file")
    subparser.add_argument("-a", "--allrecords", action='store_true',      dest="all_records", help="Process all VCF records assuming all records are snp records.")
    subparser.add_argument("-p", "--pass",       action='store_true',      dest="passfilter",  help="Process only the VCF samples or records having PASS FT element or PASS filter.  The filter element is always ignored when samples have the FT element regardless of this option.")
    subparser.add_argument("-t", "--tableFile",  type=str, metavar='FILE', dest="table_file",  help="Tablulate the results in the specified tab-separated-value file.")
    subparser.set_defaults(func=compare)

    args = parser.parse_args(system_args)
    return args


def compare(args):
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
    args.vcf_path_list = args.vcf_path1 + args.vcf_path_list

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
    snp_set_list = [set(sample_snps.snp_list) for sample_snps in sample_snps_list] # list of set of SnpTuple namedtuples having (chrom, pos, ref, gt, sample)
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
        positions = set([(s.chrom, s.pos) for s in snp_set_list[i]])
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
            print("CHROM   \tPOS")
            for chrom, pos in sorted_positions:
                print(chrom + '\t' + str(pos))

    # Print the positions missing in each of the VCF files
    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            print("\nPositions missing in %s, but present in at least 2 other VCF files:" % base_vcf_file_name_list[i])
            sorted_positions = sorted(list(missing_positions_sets[i]))
            if len(sorted_positions) == 0:
                print("None")
            else:
                print("CHROM   \tPOS")
                for chrom, pos in sorted_positions:
                    print(chrom + '\t' + str(pos))

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

    # Tabulate all the venn diagram sections into a TSV spreadsheet
    # Emit one row per position and one column per sample.
    if args.table_file:
        tabulate_results(snp_set_list, all_samples, args.table_file)

    # Generate a venn diagram if the necessary packages are installed
    #generate_venn_diagrams(snp_set_list, base_vcf_file_name_list, 'snps.venn.pdf')
    if num_vcf_files >= 2 and num_vcf_files <= 6:
        try:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib import pyplot as plt
        except:
            report_error("Skipping venn diagram creation.  Requires matplotlib.")
            exit(1)
    if num_vcf_files == 2 or num_vcf_files == 3:
        try:
            from matplotlib_venn import venn2, venn3, venn3_unweighted, venn2_circles, venn3_circles
        except:
            report_error("Skipping venn diagram creation.  Requires matplotlib_venn.")
            exit(1)

        def colorize_venn2(venn_circles):
            colors2 = ['red', 'blue', 'magenta']
            alpha2 = [0.6, 0.4, 0.1]
            venn_circles.get_patch_by_id('10').set_color(colors2[0])
            venn_circles.get_patch_by_id('01').set_color(colors2[1])
            venn_circles.get_patch_by_id('11').set_color(colors2[2])
            venn_circles.get_patch_by_id('10').set_alpha(alpha2[0])
            venn_circles.get_patch_by_id('01').set_alpha(alpha2[1])
            venn_circles.get_patch_by_id('11').set_alpha(alpha2[2])

        def make_venn2(set_list, title, output_file):
            """
            Draw multiple stacked 2-circle venn diagrams with all pairs of sets.
            """
            num_sets = len(set_list)
            if num_sets == 2:
                c = venn2(set_list, set_labels=base_vcf_file_name_list)
                colorize_venn2(c)
                plt.title(title)
            else:
                fig, axes = plt.subplots(3)
                plt_idx = 0
                for pair in itertools.combinations(range(num_sets), 2):
                    sets = [set_list[k] for k in range(num_sets) if k in pair]
                    names = [base_vcf_file_name_list[k] for k in range(num_sets) if k in pair]
                    c = venn2(sets, set_labels=names, ax=axes[plt_idx])
                    colorize_venn2(c)
                    plt_idx += 1
                axes[0].set_title(title)
            plt.show()
            plt.savefig(output_file)
            plt.close()

        make_venn2(snp_set_list, "Venn Diagram of SNPs", "venn2.snps.pdf")
        make_venn2(position_set_list, "Venn Diagram of Positions", "venn2.positions.pdf")

    if num_vcf_files >= 3 and num_vcf_files <= 6:
        from pyvenn import venn

        name_colors = ["black" for i in range(num_vcf_files)]
        figsize = (5, 8)
        fontsize = 12
        titlesize = 14
        legend_loc = None
        venn_func_dict = {3:venn.venn3, 4:venn.venn4, 5:venn.venn5, 6:venn.venn6}
        fig, axes = plt.subplots(1, 2, figsize=(11, 8), tight_layout=True)

        plt_idx = 0
        labels = venn.get_labels(position_set_list, fill=["number"])
        venn_func_dict[num_vcf_files](labels, base_vcf_file_name_list, figsize=figsize, fontsize=fontsize, name_colors=name_colors, legend_loc=legend_loc, axes=axes[plt_idx])
        axes[plt_idx].set_title("Positions", fontsize=titlesize)

        plt_idx += 1
        labels = venn.get_labels(snp_set_list, fill=["number"])
        venn_func_dict[num_vcf_files](labels, base_vcf_file_name_list, figsize=figsize, fontsize=fontsize, name_colors=name_colors, legend_loc=legend_loc, axes=axes[plt_idx])
        axes[plt_idx].set_title("SNPs", fontsize=titlesize)

        plt.subplots_adjust(hspace=0.5)  # adjust height spacing between subplots
        plt.show()
        plt.savefig("venn%i.pdf" % num_vcf_files)
        plt.close()


if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    args.func(args)

