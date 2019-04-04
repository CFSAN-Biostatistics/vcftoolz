# -*- coding: utf-8 -*-

"""This module is part of VCF Toolz.
"""

from __future__ import print_function
from __future__ import absolute_import

from Bio import SeqIO
import csv
import collections
import itertools
import logging
import numpy as np
import os
import pandas as pd
import sys
from timeit import default_timer as timer
import vcf
from vcftoolz.vcf_call_parser import call_alleles, call_generator, is_snp_call, is_indel_call, is_other_variant_call, is_ref_call, is_filtered_call, is_missing_call


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

    Examples
    --------
    get_unique_set_elements([{"100","110","101","111"}, {"110","010","111","011"}, {"101","111","011","001"}])
    [{'100'}, {'010'}, {'001'}]
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

    Examples
    --------
    get_missing_set_elements([{"100","110","101","111"}, {"110","010","111","011"}, {"101","111","011","001"}])
    [{'011'}, {'101'}, {'110'}]
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

    Examples
    --------
    get_set_intersections([{"100","110","101","111"}, {"110","010","111","011"}, {"101","111","011","001"}])
    [('111', {'111'}), ('011', {'011'}), ('101', {'101'}), ('001', {'001'}), ('110', {'110'}), ('010', {'010'}), ('100', {'100'})]
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
            logging.error(err_message)
            continue
        if os.path.getsize(file_path) == 0:
            bad_count += 1
            err_message = "%s %s is empty." % (error_prefix, file_path)
            logging.error(err_message)
            continue

    return bad_count


def open_vcf_file(vcf_path):
    """Open a vcf file for reading. Detects and handles gz compressed files.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file

    Returns
    -------
    open file handle
    """
    bad_files_count = verify_non_empty_input_files("VCF file", [vcf_path])
    if bad_files_count > 0:
        exit(1)
    compressed = vcf_path.lower().endswith(".gz")
    mode = "rb" if compressed else "rt"
    handle = open(vcf_path, mode)
    return handle


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
    reader = vcf.VCFReader(filename=vcf_path)
    return reader.samples


SnpTuple = collections.namedtuple("SnpTuple", "chrom pos ref gt sample")


class SampleSnps(object):
    """
    Class to hold a list of snps and the genotype
    data per snp position.
    """

    def __init__(self, snp_list, filtered_snp_list, alt_dict, format_dict, call_dict):
        """
        Initialize the Snps container.

        Parameters
        ----------
        snp_list : list of SnpTuple namedtuple
            List of SnpTuple namedtuple containing (chrom, pos, ref, gt, sample)
        filtered_snp_list : list of SnpTuple namedtuple
            List of filtered SnpTuple namedtuple containing (chrom, pos, ref, gt, sample)
        alt_dict :  dict
            Dictionary mapping from keyed snp tuple to actual alt list
        format_dict : dict
            Dictionary mapping from keyed snp tuple to FORMAT string
        call_dict : dict
            Dictionary mapping from keyed snp tuple to namedtuple of genotype data elements
        """
        self.snp_list = snp_list
        self.filtered_snp_list = filtered_snp_list
        self.alt_dict = alt_dict
        self.format_dict = format_dict
        self.call_dict = call_dict


def get_snp_list(truth_flag, vcf_path, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
    """
    Get the list of snps in a VCF file.

    By default, all calls are included in the output, unless exluded.

    Parameters
    ----------
    truth_flag : bool
        When true, the first VCF file is considered truth and the filtered snps are gathered separately from non-filtered snps
    vcf_path : str
        Path to the VCF file
    exclude_snps : bool
        Exclude snp calls.
    exclude_indels : bool
        Exclude insertions and deletions.
    exclude_vars : bool
        Exclude variants other than snps and indels.
    exclude_refs : bool
        Exclude reference calls.
    exclude_hetero : bool
        Exclude heterozygous calls.
    exclude_filtered : bool
        Exclude filtered calls (FT or FILTER is not PASS).
    exclude_missing : bool
        Exclude calls with all data elements missing.

    Returns
    -------
    snps : SampleSnps
        Container of:
            snps : list of namedtuple SnpTuple containing (chrom, pos, ref, gt, sample)
            filtered_snps : list of namedtuple SnpTuple containing (chrom, pos, ref, gt, sample)
            alt_dict :  dictionary mapping from keyed snp tuple to actual alt list
            format_dict : dictionary mapping from keyed snp tuple to FORMAT string
            call_dict : dictionary mapping from keyed snp tuple to namedtuple of genotype data elements
    """
    compressed = vcf_path.lower().endswith(".gz")
    mode = "rb" if compressed else "rt"
    with open(vcf_path, mode) as input:
        snps = []
        filtered_snps = []
        alt_dict = dict()
        format_dict = dict()
        call_dict = dict()

        if truth_flag:
            exclude_filtered = False

        for record, call in call_generator(input, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
            if call.is_het:
                bases = call.gt_bases or '.'  # heterozygous
            else:
                bases = call_alleles(call)[0]  # homozygous
            if bases == "N":
                continue
            snp = SnpTuple(record.CHROM, int(record.POS), record.REF.upper(), bases, call.sample)
            filtered = is_filtered_call(record, call)
            if not truth_flag or not filtered:
                snps.append(snp)
            if filtered:
                filtered_snps.append(snp)
            alt_dict[snp] = record.ALT
            format_dict[snp] = record.FORMAT
            call_dict[snp] = call.data
        return SampleSnps(snps, filtered_snps, alt_dict, format_dict, call_dict)


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
            positions = sorted(positions)  # set becomes sorted list
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


def is_list_like(obj):
    """Determine if an object is an iterable collection, but not a string.

    Parameters
    ----------
    obj : object
        Object to inspect.

    Returns
    -------
    true if object is a collection, but not a string.

    Examples
    --------
    >>> is_list_like(None)
    False
    >>> is_list_like(0)
    False
    >>> is_list_like("aaa")
    False
    >>> is_list_like(r"aaa")
    False
    >>> is_list_like(u"aaa")
    False
    >>> is_list_like(b"aaa")
    False
    >>> is_list_like([1,2,3])
    True
    >>> is_list_like((1,2,3))
    True
    """
    if obj is None:
        return False
    if isinstance(obj, str):
        return False
    if isinstance(obj, bytes):
        return False
    if hasattr(obj, "__iter__"):
        return True
    return False


def stringify_call_data_list(call_data_list, separator):
    """Convert None to '.' while joining items in a list to make a string.

    Parameters
    ----------
    call_data_list : list
        VCF call data list
    separator : str
        Separator for items in list.

    Returns
    -------
    s : str
        List items joined with None converted to '.'
    """
    call_data_list = ['.' if v is None else str(v) for v in call_data_list]
    return separator.join(call_data_list)


def print_snp_detail_list(snp_set, alt_dict, format_dict, call_dict):
    """Print a list of snps.

    Parameters
    ---------
    snp_set : list of SnpTuple namedtuple
        Set of SnpTuple namedtuple containing (chrom, pos, ref, gt, sample)
    alt_dict :  dict
        Dictionary mapping from keyed snp tuple to actual alt list
    format_dict : dict
        Dictionary mapping from keyed snp tuple to FORMAT string
    call_dict : dict
        Dictionary mapping from keyed snp tuple to namedtuple of genotype data elements
    """
    sorted_snps = sorted(list(snp_set))

    if len(sorted_snps) == 0:
        print("None")
    else:
        print("CHROM   \tPOS\tREF\tALT\tSAMPLE  \tFORMAT\tDATA")
        for snp in sorted_snps:
            fields = [str(x) for x in snp]
            fields[3] = ",".join(map(str, alt_dict[snp]))
            format_str = format_dict[snp]
            format_keys = format_str.split(":")
            call_data = call_dict[snp]
            call_data_list = [getattr(call_data, k, None) for k in format_keys]
            call_data_list = [stringify_call_data_list(v, ',') if is_list_like(v) else v for v in call_data_list]
            call_data_str = stringify_call_data_list(call_data_list, ':')
            fields.append(format_str)
            fields.append(call_data_str)
            print('\t'.join(fields))


def compare(truth_flag, vcf_path_list, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing, table_file):
    """
    Compare and analyze the snps found in two or more input VCF files.

    By default, all calls are included in the output.

    Parameters
    ----------
    truth_flag : bool
        When true, the first VCF file is considered truth.
    vcf_path_list : list of str
        List of paths to the VCF files
    exclude_snps : bool
        Exclude snp calls.
    exclude_indels : bool
        Exclude insertions and deletions.
    exclude_vars : bool
        Exclude variants other than snps and indels.
    exclude_refs : bool
        Exclude reference calls.
    exclude_hetero : bool
        Exclude heterozygous calls.
    exclude_filtered : bool
        Exclude filtered calls (FT or FILTER is not PASS).
    exclude_missing : bool
        Exclude calls with all data elements missing.
    table_file_path : str
        Path to the TSV file to be written.
    """
    # Validate input files
    bad_files_count = verify_non_empty_input_files("VCF file", vcf_path_list)
    if bad_files_count > 0:
        exit(1)

    # Get short base file name for each input
    base_vcf_file_name_list = [os.path.basename(p) for p in vcf_path_list]
    base_vcf_file_name_list = [f[: -3] if f.lower().endswith(".gz") else f for f in base_vcf_file_name_list]
    base_vcf_file_name_list = [f[: -4] if f.lower().endswith(".vcf") else f for f in base_vcf_file_name_list]
    num_vcf_files = len(base_vcf_file_name_list)

    # Get the list of sample names in each VCF file
    sample_set_list = [set(get_sample_list(path)) for path in vcf_path_list]
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
    packed_exclude_flags = (exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing)
    sample_snps_list = [get_snp_list(truth_flag, path, *packed_exclude_flags) for path in vcf_path_list]  # List of SampleSnps
    snp_set_list = [set(sample_snps.snp_list) for sample_snps in sample_snps_list]  # list of set of SnpTuple namedtuples having (chrom, pos, ref, gt, sample)
    filtered_snp_set_list = [set(sample_snps.filtered_snp_list) for sample_snps in sample_snps_list]  # list of set of SnpTuple namedtuples having (chrom, pos, ref, gt, sample)
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
    filtered_position_set_list = []
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        positions = set([(s.chrom, s.pos) for s in snp_set_list[i]])
        filtered_positions = set([(s.chrom, s.pos) for s in filtered_snp_set_list[i]])
        position_set_list.append(positions)
        filtered_position_set_list.append(filtered_positions)
        count = len(positions)
        print("{count}\tPositions having snps in {dataset}".format(dataset=dataset, count=count))
    print()

    unique_positions_sets = get_unique_set_elements(position_set_list)
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(unique_positions_sets[i])
        if not truth_flag:
            unique_snps_description = "SNP positions only in {dataset}"
        elif i == 0:
            unique_snps_description = "False negative snp positions"
        else:
            unique_snps_description = "False positive snp positions"
        print(("{count}\t" + unique_snps_description).format(dataset=dataset, count=count))
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
        print("{count}\tCalls in {dataset}".format(dataset=dataset, count=count))
    print()

    if truth_flag:
        true_positive_snps = snp_set_list[0] & snp_set_list[1]
        true_positive_snps_count = len(true_positive_snps)
        print("{count}\tTrue positive calls".format(count=true_positive_snps_count))
    for i in range(num_vcf_files):
        dataset = base_vcf_file_name_list[i]
        count = len(unique_snps_sets[i])
        if not truth_flag:
            unique_snps_description = "Calls only in {dataset}"
        elif i == 0:
            unique_snps_description = "False negative calls"
            false_negative_snps_count = count
        else:
            unique_snps_description = "False positive calls"
            false_positive_snps_count = count
        print(("{count}\t" + unique_snps_description).format(dataset=dataset, count=count))
    if truth_flag:
        precision = float(true_positive_snps_count) / (true_positive_snps_count + false_positive_snps_count)
        recall = float(true_positive_snps_count) / (true_positive_snps_count + false_negative_snps_count)
        f1 = 2.0 * (precision * recall) / (precision + recall)
        print("%0.4f\tPrecision" % precision)
        print("%0.4f\tRecall" % recall)
        print("%0.4f\tF1 score" % f1)
    print()

    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            dataset = base_vcf_file_name_list[i]
            count = len(missing_snps_sets[i])
            print("{count}\tCalls missing in {dataset}, but present in at least 2 other VCF files".format(dataset=dataset, count=count))
        print()

    # Print the effectiveness of the snp filters
    if truth_flag:
        truth_positions = position_set_list[0]
        filtered_positions = filtered_position_set_list[1]
        incorrectly_filtered_positions = truth_positions & filtered_positions
        correctly_filtered_positions = filtered_positions - truth_positions

        # Need to do something special here and omit the gt element when examining filtered snps
        # because the gt element is usually "." which will never match the truth gt element. So
        # instead, we just look at chrom, pos, and sample.
        truth_snps = set([(s.chrom, s.pos, s.sample) for s in snp_set_list[0]])
        filtered_snps = set([(s.chrom, s.pos, s.sample) for s in filtered_snp_set_list[1]])
        incorrectly_filtered_snps = truth_snps & filtered_snps
        correctly_filtered_snps = filtered_snps - truth_snps

        print("{count}\tTotal filtered snp positions".format(count=len(filtered_positions)))
        print("{count}\tFalse negative snp positions incorrectly removed by filters".format(count=len(incorrectly_filtered_positions)))
        print("{count}\tTrue negative snp positions correctly removed by filters".format(count=len(correctly_filtered_positions)))
        print()
        print("{count}\tTotal filtered calls".format(count=len(filtered_snps)))
        print("{count}\tFalse negative calls incorrectly removed by filters".format(count=len(incorrectly_filtered_snps)))
        print("{count}\tTrue negative calls correctly removed by filters".format(count=len(correctly_filtered_snps)))
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
        print("\nCalls only in %s" % base_vcf_file_name_list[i])
        print_snp_detail_list(unique_snps_sets[i], alt_dict_list[i], format_dict_list[i], call_dict_list[i])

    # Print the snps missing in each of the VCF files
    if num_vcf_files >= 3:
        for i in range(num_vcf_files):
            print("\nCalls missing in %s, but present in at least 2 other VCF files:" % base_vcf_file_name_list[i])
            sorted_snps = sorted(list(missing_snps_sets[i]))
            if len(sorted_snps) == 0:
                print("None")
            else:
                print("CHROM   \tPOS\tREF\tALT\tSAMPLE")
                for snp in sorted_snps:
                    fields = [str(x) for x in snp]
                    print('\t'.join(fields))

    # Print the false negative filtered snps
    if truth_flag:
        # Get the full SnpTuple for all known incorrectly filtered snps
        incorrectly_filtered_snps = set([s for s in filtered_snp_set_list[1] if (s.chrom, s.pos, s.sample) in incorrectly_filtered_snps])
        print("\nFalse negative calls incorrectly removed by filters")
        print_snp_detail_list(incorrectly_filtered_snps, alt_dict_list[1], format_dict_list[1], call_dict_list[1])

    # Print the true negative filtered snps
    if truth_flag:
        # Get the full SnpTuple for all known correctly filtered snps
        correctly_filtered_snps = set([s for s in filtered_snp_set_list[1] if (s.chrom, s.pos, s.sample) in correctly_filtered_snps])
        print("\nTrue negative calls correctly removed by filters")
        print_snp_detail_list(correctly_filtered_snps, alt_dict_list[1], format_dict_list[1], call_dict_list[1])

    # Tabulate all the venn diagram sections into a TSV spreadsheet
    # Emit one row per position and one column per sample.
    if table_file:
        tabulate_results(snp_set_list, all_samples, table_file)

    # Generate a venn diagram if the necessary packages are installed
    if num_vcf_files >= 2 and num_vcf_files <= 6:
        try:
            import matplotlib
            matplotlib.use("Agg")
            from matplotlib import pyplot as plt
        except ImportError:
            logging.info("Skipping venn diagram creation.  Requires matplotlib.")
            exit(1)
    if num_vcf_files == 2 or num_vcf_files == 3:
        try:
            from matplotlib_venn import venn2, venn3, venn3_unweighted, venn2_circles, venn3_circles  # noqa F401
        except ImportError:
            logging.info("Skipping venn diagram creation.  Requires matplotlib_venn.")
            exit(1)

        def colorize_venn2(venn_circles):
            """Set the colors of the 3 sections of a 2-set Venn diagram.
            """
            colors2 = [
                # r, g, b
                [157, 217, 161],
                [156, 195, 229],
                [117, 180, 192],
            ]
            colors2 = [[i[0] / 255.0, i[1] / 255.0, i[2] / 255.0] for i in colors2]

            for idx, patch_id in enumerate(["10", "01", "11"]):
                patch = venn_circles.get_patch_by_id(patch_id)
                if patch:
                    patch.set_color(colors2[idx])
                    patch.set_alpha(1.0)

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
        venn_func_dict = {3: venn.venn3, 4: venn.venn4, 5: venn.venn5, 6: venn.venn6}
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
        plt.savefig("venn%i.pdf" % num_vcf_files)
        plt.close()


def build_narrow_row(record, call):
    """Construct a list of narrow format row items for a single call in the narrow format.

    Parameters
    ----------
    record : pyvcf Record
        Example: Record(CHROM=CP006053.1, POS=3636502, REF=T, ALT=[N, G])
    call : pyvcf Call
        Example: Call(sample=CFSAN000212, CallData(GT=2, SDP=8, RD=0, AD=[3, 5], RDF=0, RDR=0, ADF=[1, 5], ADR=[2, 0], FT=PASS))
        Example: Call(sample=CFSAN000191, CallData(GT=., SDP=2, RD=1, AD=[1, None], RDF=0, RDR=1, ADF=[1, None], ADR=[0, None], FT=VarFreq60;Depth3))

    Returns
    -------
    row : list
        List of narrow format row items for a single call.

    Examples
    --------
    >>> call_format = "GT:SDP:RD:AD:RDF:RDR:ADF:ADR:FT"
    >>> CallData = vcf.model.make_calldata_tuple(call_format.split(':'))

    >>> # Single alternate allele, GT=1
    >>> # CallData(GT=1, SDP=3, RD=0, AD=[3], RDF=0, RDR=0, ADF=[1], ADR=[2], FT=PASS))
    >>> record = vcf.model._Record("CP006053.1", 3636502, "ID", "T", [vcf.model._Substitution(b) for b in ['G']], '.', "PASS", "NS=1", call_format, None)
    >>> sample_data = CallData('1', 3, 0, [3], 0, 0, [1], [2], "PASS")
    >>> call = vcf.model._Call(record, "CFSAN000191", sample_data)
    >>> build_narrow_row(record, call)
    ['CFSAN000191', 'CP006053.1', 3636502, 'T', 'G', '1', 3, 0, '3', 0, 0, '1', '2', 'PASS']

    >>> # Single alternate allele, GT=.
    >>> # CallData(GT=., SDP=3, RD=0, AD=[3], RDF=0, RDR=0, ADF=[1], ADR=[2], FT=PASS))
    >>> record = vcf.model._Record("CP006053.1", 3636502, "ID", "T", [vcf.model._Substitution(b) for b in ['G']], '.', "PASS", "NS=1", call_format, None)
    >>> sample_data = CallData('.', 3, 0, [3], 0, 0, [1], [2], "PASS")
    >>> call = vcf.model._Call(record, "CFSAN000191", sample_data)
    >>> build_narrow_row(record, call)
    ['CFSAN000191', 'CP006053.1', 3636502, 'T', 'G', '.', 3, 0, '3', 0, 0, '1', '2', 'PASS']

    >>> # Multiple alternate alleles, GT=2
    >>> # CallData(GT=2, SDP=8, RD=0, AD=[3, 5], RDF=0, RDR=0, ADF=[1, 5], ADR=[2, 0], FT=PASS))
    >>> record = vcf.model._Record("CP006053.1", 3636502, "ID", "T", [vcf.model._Substitution(b) for b in ['N', 'G']], '.', "PASS", "NS=1", call_format, None)
    >>> sample_data = CallData('2', 8, 0, [3, 5], 0, 0, [1, 5], [2, 0], "PASS")
    >>> call = vcf.model._Call(record, "CFSAN000191", sample_data)
    >>> build_narrow_row(record, call)
    ['CFSAN000191', 'CP006053.1', 3636502, 'T', 'N,G', '2', 8, 0, '3,5', 0, 0, '1,5', '2,0', 'PASS']

    >>> # Multiple alternate alleles, GT=.
    >>> # CallData(GT=., SDP=2, RD=1, AD=[1, None], RDF=0, RDR=1, ADF=[1, None], ADR=[0, None], FT=VarFreq60;Depth3))
    >>> record = vcf.model._Record("CP006053.1", 3636502, "ID", "T", [vcf.model._Substitution(b) for b in ['N', 'G']], '.', "PASS", "NS=1", call_format, None)
    >>> sample_data = CallData('.', 2, 1, [1, None], 0, 1, [1, None], [0, None], "VarFreq60;Depth3")
    >>> call = vcf.model._Call(record, "CFSAN000191", sample_data)
    >>> build_narrow_row(record, call)
    ['CFSAN000191', 'CP006053.1', 3636502, 'T', 'N,G', '.', 2, 1, '1,.', 0, 1, '1,.', '0,.', 'VarFreq60;Depth3']
    """
    alt_list = [str(variant) for variant in record.ALT]
    alt_bases_str = ','.join(alt_list)
    call_data_list = [stringify_call_data_list(v, ',') if is_list_like(v) else v for v in call.data]
    call_data_list = ['.' if item is None else item for item in call_data_list]
    row = [call.sample, record.CHROM, int(record.POS), record.REF, alt_bases_str] + call_data_list
    return row


def narrow(vcf_path, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
    """Convert a VCF file into a tab delimited set of snp calls, one per line.

    By default, all calls are included in the output.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    exclude_snps : bool
        Exclude snp calls.
    exclude_indels : bool
        Exclude insertions and deletions.
    exclude_vars : bool
        Exclude variants other than snps and indels.
    exclude_refs : bool
        Exclude reference calls.
    exclude_hetero : bool
        Exclude heterozygous calls.
    exclude_filtered : bool
        Exclude filtered calls (FT or FILTER is not PASS).
    exclude_missing : bool
        Exclude calls with all data elements missing.
    """
    input = open_vcf_file(vcf_path)

    snps = []
    for record, call in call_generator(input, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
        row = build_narrow_row(record, call)
        snps.append(row)

    snps = sorted(snps, key=lambda snp: snp[1:3] + snp[0:1])  # sort by chrom, then pos, then sample

    out = csv.writer(sys.stdout, delimiter='\t', lineterminator=os.linesep)
    header = ["Sample", "CHROM", "POS", "REF", "ALT"]
    if record:
        header += record.FORMAT.split(':')
    out.writerow(header)
    for row in snps:
        out.writerow(row)


def count(vcf_path, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
    """Count the number of positions and calls.

    By default, all calls are included in the output.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    exclude_snps : bool
        Exclude snp calls.
    exclude_indels : bool
        Exclude insertions and deletions.
    exclude_vars : bool
        Exclude variants other than snps and indels.
    exclude_refs : bool
        Exclude reference calls.
    exclude_hetero : bool
        Exclude heterozygous calls.
    exclude_filtered : bool
        Exclude filtered calls (FT or FILTER is not PASS).
    exclude_missing : bool
        Exclude calls with all data elements missing.
    """
    input = open_vcf_file(vcf_path)

    position_set = set()
    call_count = 0
    heterozygous_count = 0
    homozygous_count = 0
    ref_count = 0
    snp_count = 0
    indel_count = 0
    other_variant_count = 0
    filtered_count = 0
    passing_count = 0
    missing_count = 0
    filtered_reason_counts = collections.Counter()
    for record, call in call_generator(input, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
        position_set.add(record.POS)
        call_count += 1

        if call.is_het:
            heterozygous_count += 1
        else:
            homozygous_count += 1

        if is_ref_call(record, call):
            ref_count += 1
        if is_snp_call(record, call):
            snp_count += 1
        if is_indel_call(record, call):
            indel_count += 1
        if is_other_variant_call(record, call):
            other_variant_count += 1
        if is_missing_call(record, call):
            missing_count += 1

        try:
            filt = call.data.FT     # If there are filters per sample, use them
        except AttributeError:
            filt = record.FILTER    # Otherwise, use the filter for the whole record
        if filt is None:
            pass
        elif filt == "PASS":
            passing_count += 1
        else:
            filtered_count += 1

            # There might be multiple filter reasons per call
            filters = filt.split(';')
            for filter in filters:
                filtered_reason_counts[filter] += 1

    sample_list = get_sample_list(vcf_path)

    rows = []
    rows.append(("samples", len(sample_list)))
    rows.append(("positions", len(position_set)))
    rows.append(("calls", call_count))
    rows.append(("heterozygous calls", heterozygous_count))
    rows.append(("homozygous calls", homozygous_count))
    rows.append(("reference calls", ref_count))
    rows.append(("snp calls", snp_count))
    rows.append(("indel calls", indel_count))
    rows.append(("other variant calls", other_variant_count))
    rows.append(("missing calls", missing_count))
    rows.append(("calls PASS filter", passing_count))
    rows.append(("calls failing any filter", filtered_count))
    for filter, count in filtered_reason_counts.most_common():
        rows.append(("calls failing filter %s" % filter, count))
    max_width = 0
    for label, count in rows:
        max_width = max(max_width, len(str(count)))
    for label, count in rows:
        print("%*d %s" % (max_width, count, label))


def plot(vcf_path, reference_path, output_path, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing, chrom=None):
    """Count the number of positions and calls.

    By default, all calls are included in the output.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file
    reference_path : str
        Path to the reference fasta file
    output_path : str
        Output file.  You should use the extension .pdf or .png.
    exclude_snps : bool
        Exclude snp calls.
    exclude_indels : bool
        Exclude insertions and deletions.
    exclude_vars : bool
        Exclude variants other than snps and indels.
    exclude_refs : bool
        Exclude reference calls.
    exclude_hetero : bool
        Exclude heterozygous calls.
    exclude_filtered : bool
        Exclude filtered calls (FT or FILTER is not PASS).
    exclude_missing : bool
        Exclude calls with all data elements missing.
    chrom : str, optional, defaults to None
        When specified, restricts the plot to a single contig.  Otherwise, all the contigs are joined together
        in order of decreasing size.
    """
    import matplotlib.pylab as plt

    # Determine the length of each contig
    contig_length_dict = collections.Counter()
    with open(reference_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if chrom is None or chrom == record.id:
                contig_length_dict[record.id] = len(record.seq)

    # Collect all the records and calls. Store them so we can loop through the list more than once below.
    start = timer()
    input = open_vcf_file(vcf_path)
    record_call_tuples = list()
    for record, call in call_generator(input, exclude_snps, exclude_indels, exclude_vars, exclude_refs, exclude_hetero, exclude_filtered, exclude_missing):
        if chrom is None or chrom == record.CHROM:
            record_call_tuples.append((record, call))
    input.close()
    end = timer()
    logging.info("%.1f seconds parsing vcf file" % (end - start))

    # Count calls passing all filters
    # defaultdict key = contig, value = Counter of calls per position
    start = timer()
    contig_call_counts = collections.defaultdict(collections.Counter)
    for record, call in record_call_tuples:
        if not is_filtered_call(record, call):
            contig_call_counts[record.CHROM][record.POS] += 1
    end = timer()
    logging.info("%.1f seconds counting calls passing all filters" % (end - start))

    # Count calls failing filters
    # defaultdict key = filter, value = (defaultdict key = contig, value = Counter of calls per position)
    start = timer()
    filtered_reason_counts = collections.Counter()
    filtered_contig_call_counts = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for record, call in record_call_tuples:
        try:
            filt = call.data.FT     # If there are filters per sample, use them
        except AttributeError:
            filt = record.FILTER    # Otherwise, use the filter for the whole record
        if not filt or filt == "PASS":  # filter is not set or set to PASS
            continue

        # There might be multiple filter reasons per call
        filters = filt.split(';')
        for filter in filters:
            filtered_reason_counts[filter] += 1
            filtered_contig_call_counts[filter][record.CHROM][record.POS] += 1
    end = timer()
    logging.info("%.1f seconds counting calls failing filters" % (end - start))

    def plotter(ax_left, ax_right, contig_length_dict, contig_call_counts, title):
        # Create an arbitrary genome with the contigs lined-up in order of decreasing size.
        # Count the calls at the genome positions.
        length_of_all_prior_contigs = 0
        genome_calls_per_pos = dict()
        for contig, length in contig_length_dict.most_common():
            calls_per_contig_pos = contig_call_counts[contig]
            for contig_pos in calls_per_contig_pos:
                genome_calls_per_pos[contig_pos + length_of_all_prior_contigs] = calls_per_contig_pos[contig_pos]
            length_of_all_prior_contigs += length

        # Put entries at first and last position to make sure cumsum plot extends all the way left and right
        if 1 not in genome_calls_per_pos:
            genome_calls_per_pos[1] = 0
        if length_of_all_prior_contigs not in genome_calls_per_pos:
            genome_calls_per_pos[length_of_all_prior_contigs] = 0
        lists = sorted(genome_calls_per_pos.items())  # list of tuples, sorted by position
        x, y = zip(*lists)  # unpack a list of pairs into two tuples

        y_right = np.cumsum(y)
        ax_right.plot(x, y_right, color="C1", zorder=1)
        ax_right.set_ylabel("Cumulative Calls", color="C1")

        ax_left.bar(x, y, linewidth=1, edgecolor="C0", zorder=2)
        ax_left.set_ylabel("Calls", color="C0")
        ax_left.set_title(title)

        return np.max(x), np.max(y), np.max(y_right)

    start = timer()
    num_subplots = 1 + len(filtered_reason_counts)
    figure, axes_left = plt.subplots(nrows=num_subplots, ncols=1, figsize=(8.5, 2 * num_subplots), tight_layout=True)
    if num_subplots == 1:
        axes_left = [axes_left]
    axes_right = [ax.twinx() for ax in axes_left]
    title = "Calls Passing All Filters"
    xmax, ymax_left, ymax_right = plotter(axes_left[0], axes_right[0], contig_length_dict, contig_call_counts, title)
    call_counts = [ymax_right]
    call_reasons = [title]
    for i, (filter, _) in enumerate(filtered_reason_counts.most_common()):
        contig_call_counts = filtered_contig_call_counts[filter]
        title = "Calls Failing Filter " + filter
        ax_xmax, ax_left_ymax, ax_right_ymax = plotter(axes_left[1 + i], axes_right[1 + i], contig_length_dict, contig_call_counts, title)
        xmax = max(xmax, ax_xmax)
        ymax_left = max(ymax_left, ax_left_ymax)
        ymax_right = max(ymax_right, ax_right_ymax)
        call_counts.append(ax_right_ymax)
        call_reasons.append(title)
    xmax *= 1.00
    ymax_left *= 1.07
    ymax_right *= 1.07
    for ax in axes_left:
        ax.set_xlim(left=0, right=xmax, auto=True)
        ax.set_ylim(bottom=0, top=ymax_left, auto=True)
    for ax in axes_right:
        ax.set_ylim(bottom=0, top=ymax_right, auto=True)

    plt.savefig(output_path)
    plt.close()
    end = timer()
    logging.info("%.1f seconds plotting" % (end - start))

    # Print table of contig start-end positions
    df = pd.DataFrame(contig_length_dict.most_common())
    df.columns = ["chrom", "length"]
    df["end"] = np.cumsum(df.length)
    df["start"] = 1 + df.end - df.length
    df = df[["chrom", "start", "end", "length"]]
    df.columns = ["Chrom", "Start", "End", "Length"]
    print(df.to_string(index=False))
    print()

    # Print table of call counts
    df = pd.DataFrame()
    df["Filter"] = call_reasons
    df["Count"] = call_counts
    print(df.to_string(index=False))
