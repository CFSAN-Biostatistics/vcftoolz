# -*- coding: utf-8 -*-

"""
Module to parse vcf files and generate calls with their containing records.
The records and calls are in PyVCF format.

How to test this module:
    python -m doctest vcf_call_parser.py
"""

from __future__ import print_function

import sys
import vcf


def _make_test_vcf_file(test_record):
    """This is a helper function used to help test various types of calls with the doctest module.

    It creates an in-memory VCF file containing the specified VCF record and returns a file
    handle ready for parsing by PyVCF.

    Parameters
    ----------
    test_record : str
        VCF record just as it would appear in a VCF file.
        It may contain multiple genoptype calls.

    Returns
    -------
    Open file handle positioned at the start of the in-memory VCF file containing a single
    VCF record.
    """
    if sys.version_info < (3,):
        from StringIO import StringIO
    else:
        from io import StringIO
    in_memory_file = StringIO()
    in_memory_file.name = "test.vcf"
    fields = test_record.split()
    num_fields = len(fields)
    if num_fields < 10:
        raise ValueError("Test VCF record has only %i fields, expecting at least 10" % num_fields)
    if num_fields == 10:
        sample_names_str = "SAMPLE"
    else:
        sample_names_str = ' '.join("SAMPLE%i" % i for i in range(1, 2+num_fields-10))
    header_line = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT " + sample_names_str
    print(header_line, file=in_memory_file)
    print(test_record, file=in_memory_file)
    in_memory_file.seek(0)
    return in_memory_file


def _make_test_pyvcf_calls(test_record):
    """This is a helper function used to help test various types of calls with the doctest module.

    Returns a list of PyVCF records and calls parsed from a specified test VCF record.

    Parameters
    ----------
    test_record : str
        VCF record just as it would appear in a VCF file.
        It may contain multiple genoptype calls.

    Returns
    -------
    list of (record, call) tuples parsed by PyVCF.
    """
    in_memory_file = _make_test_vcf_file(test_record)
    reader = vcf.VCFReader(in_memory_file)
    calls = []
    for record in reader:
        for call in record.samples:
            calls.append((record, call))
    return calls


def call_alleles(call):
    """Return a list of the bases in each allele of the specified call.

    Parameters
    ----------
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    alleles : list of str
        List of bases with one entry per allele.  For a diploid call, there will be a list of two strings.
        It is possible (and likely) there may be duplicate strings in the list.
    """
    if call.gt_bases is None:
        return ['.']

    alleles = call.gt_bases.split(call.gt_phase_char())
    return alleles


def is_snp_call(record, call):
    """Return True if a sample call is a snp.

    A snp is a call where REF is a single base and at least one GT allele is different from the REF

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is a snp.

    Examples
    --------
    >>> # normal snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    True
    >>> # an insertion is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A AT . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a deletion is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a homozygous spanning deletion (*) is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A * . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a heterozygous spanning deletion (*) and snp together is a snp because at least one is a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A *,G . PASS . GT:FT 1/2:PASS")[0])
    True
    >>> # a heterozygous snp and insertion together is a snp because at least one is a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A G,AT . PASS . GT:FT 1/2:PASS")[0])
    True
    >>> # a structural variant is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS")[0])
    False
    >>> # a ref call is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 0/0:PASS")[0])
    False
    >>> # a missing call is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A . . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not a snp
    >>> is_snp_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT .:.")[0])
    False
    """
    if len(record.REF) > 1:
        return False

    if is_ref_call(record, call):
        return False

    if is_missing_call(record, call):
        return False

    # At least one snp
    if any([bases in ['A', 'C', 'G', 'T', 'N'] for bases in call_alleles(call)]):
        return True

    return False


def is_indel_call(record, call):
    """Return True if a sample call is an insertion or deletion.

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is an indel.

    Examples
    --------
    >>> # a normal snp is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # insertion
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A AT . PASS . GT:FT 1/1:PASS")[0])
    True
    >>> # deletion
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . AA A . PASS . GT:FT 1/1:PASS")[0])
    True
    >>> # a homozygous spanning deletion (*) is an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A * . PASS . GT:FT 1/1:PASS")[0])
    True
    >>> # a heterozygous spanning deletion (*) and snp together is an indel because at least one is an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A *,G . PASS . GT:FT 1/2:PASS")[0])
    True
    >>> # a heterozygous snp and insertion together is an indel because at least one is an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A G,AT . PASS . GT:FT 1/2:PASS")[0])
    True
    >>> # a structural variant is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS")[0])
    False
    >>> # a ref call is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 0/0:PASS")[0])
    False
    >>> # a missing call is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A . . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not an indel
    >>> is_indel_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT .:.")[0])
    False
    """
    if is_missing_call(record, call):
        return False

    alleles = call_alleles(call)

    # At least one spanning deletion
    if any([bases == '*' for bases in alleles]):
        return True

    if not record.is_indel:
        return False

    if record.is_sv:
        return False

    if len(record.REF) > 1:
        return True

    # At least one indel
    if any([len(bases) > 1 for bases in alleles]):
        return True

    return False


def is_other_variant_call(record, call):
    """Return True if a sample call is a variant, but not a snp or indel.

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is a variant, but not a snp or indel.

    Examples
    --------
    >>> # a normal snp is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # an insertion is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A AT . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a deletion is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . AA A . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a homozygous spanning deletion (*) is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A * . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a heterozygous spanning deletion (*) and snp together is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A *,G . PASS . GT:FT 1/2:PASS")[0])
    False
    >>> # structural variant
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS")[0])
    True
    >>> # a ref call is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 0/0:PASS")[0])
    False
    >>> # a missing call is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A . . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not an "other variant"
    >>> is_other_variant_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT .:.")[0])
    False
    """
    if is_snp_call(record, call):
        return False
    if is_indel_call(record, call):
        return False
    if call.is_variant:  # at least one GT allele is different from the REF
        # is_variant can be:
        # None  : gt == .
        # True  : gt > 0
        # False : gt == 0
        return True
    return False


def is_ref_call(record, call):
    """Return True if a sample call is a reference call.

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is a reference call.

    Examples
    --------
    >>> # a normal snp is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # an insertion is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A AT . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a deletion is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . AA A . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a homozygous spanning deletion (*) is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A * . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a heterozygous spanning deletion (*) and snp together is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A *,G . PASS . GT:FT 1/2:PASS")[0])
    False
    >>> # a structural variant is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS")[0])
    False
    >>> # ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 0/0:PASS")[0])
    True
    >>> # a missing call is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A . . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT ./.:.")[0])
    False
    >>> # a missing call is not a ref call
    >>> is_ref_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT .:.")[0])
    False
    """
    return call.gt_type == 0  # homozygous ref


def is_filtered_call(record, call):
    """Return True if a sample call is filtered (not PASS).

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is filtered.

    Examples
    --------
    >>> # FT and FILTER both PASS
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # FT PASS, FILTER not PASS, FT takes priority
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . SOMEFAIL . GT:FT 1/1:PASS")[0])
    False
    >>> # FT FAIL, FILTER PASS, FT takes priority
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:FAIL")[0])
    True
    >>> # FT ., FILTER ., considered a PASS
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . . . GT:FT 1/1:.")[0])
    False
    >>> # No FT, FILTER PASS
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT 1/1")[0])
    False
    >>> # No FT, FILTER FAIL
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . FAIL . GT 1/1")[0])
    True
    >>> # No FT, FILTER ., considered a PASS
    >>> is_filtered_call(*_make_test_pyvcf_calls("chrom 1 . A T . . . GT 1/1")[0])
    False
    """
    try:
        filt = call.data.FT     # If there are filters per sample, use them
    except AttributeError:
        filt = record.FILTER    # Otherwise, use the filter for the whole record

    if not filt or filt == "PASS":  # filter is not set or set to PASS
        return False
    else:
        return True


def is_missing_call(record, call):
    """Return True if all the genotype data is missing.

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if all the call data is missing

    Examples
    --------
    >>> # a normal snp is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # an insertion is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A AT . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a deletion is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . AA A . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a homozygous spanning deletion (*) is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A * . PASS . GT:FT 1/1:PASS")[0])
    False
    >>> # a heterozygous spanning deletion (*) and snp together is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A *,G . PASS . GT:FT 1/2:PASS")[0])
    False
    >>> # a structural variant is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS")[0])
    False
    >>> # a ref call is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT 0/0:PASS")[0])
    False
    >>> # a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A . . PASS . GT:FT ./.:.")[0])
    True
    >>> # a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT ./.:.")[0])
    True
    >>> # a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:FT .:.")[0])
    True
    >>> # partial data is not a missing call
    >>> is_missing_call(*_make_test_pyvcf_calls("chrom 1 . A T . PASS . GT:DP:FT .:30:.")[0])
    False
    """
    return call.gt_type is None and not any(call.data[1:])


def _test_call_generator(test_record, exclude_snps=None, exclude_indels=None, exclude_vars=None, exclude_refs=None, exclude_hetero=None, exclude_filtered=None, exclude_missing=None):
    """This is a helper function used to help test parsing and filtering various types of calls with the doctest module.

    It creates an in-memory VCF file containing the specified VCF record and then parses the file with the specified
    exclusion flags.  The records and calls passing the exclusion filters are printed for the doctest module to inspect.
    """
    in_memory_file = _make_test_vcf_file(test_record)
    gen = call_generator(in_memory_file, exclude_snps=exclude_snps, exclude_indels=exclude_indels, exclude_vars=exclude_vars, exclude_refs=exclude_refs, exclude_hetero=exclude_hetero, exclude_filtered=exclude_filtered, exclude_missing=exclude_missing)
    for record, sample in gen:
        print(record, sample)


def call_generator(input, exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=False):
    """Parse a VCF file and return calls with the containing record for samples matching the specified criteria.

    Parameters
    ----------
    input : file handle
        Open file handle for reading
    exclude_snps : bool
        Exclude calls with single base substitutions,
        i.e. where REF is a single base and at least one GT allele is different from the REF
    exclude_indels : bool
        Exclude calls with insertions or deletions, i.e. where the length of REF and ALT are different
    exclude_vars : bool
        Exclude calls with other variants (not snp or indel)
    exclude_refs : bool
        Exclude reference calls
    exclude_hetero : bool
        Exclude heterozygous calls
    exclude_filtered : bool
        Exclude filtered calls (not PASS)
    exclude_missing : bool
        Exclude sample calls where all the genotype data is missing

    Returns
    -------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Examples
    --------
    >>> # snp when snps included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # snp when snps excluded
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # indel when snps included
    >>> _test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # ref call when snps included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # insertion when indels included
    >>> _test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[AT]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # insertion when indels excluded
    >>> _test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # deletion when indels included
    >>> _test_call_generator("chrom 1 . AT A . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=AT, ALT=[A]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # deletion when indels excluded
    >>> _test_call_generator("chrom 1 . AT A . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # structural variant when "other variants" included
    >>> _test_call_generator("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=False, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=AT, ALT=[A]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # structural variant when "other variants" excluded
    >>> _test_call_generator("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # ref call when ref calls included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=False, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=0/0, FT=PASS))
    >>> # ref call when ref calls excluded
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # heterozygous snp call when snps and heterozygous calls included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 0/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=0/1, FT=PASS))
    >>> # heterozygous snp call when snps included, but heterozygous calls excluded
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 0/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # FT filtered snp when snps and filtered calls included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:FAIL", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=False, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1, FT=FAIL))
    >>> # FT filtered snp when snps included, but filtered calls excluded
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:FAIL", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # FILTER filtered snp when snps and filtered calls included
    >>> _test_call_generator("chrom 1 . A T . FAIL . GT 1/1", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=False, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1))
    >>> # FILTER filtered snp when snps included, but filtered calls excluded
    >>> _test_call_generator("chrom 1 . A T . FAIL . GT 1/1", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # missing call when missing calls included
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:DP:FT ./.:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=False)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=./., DP=None, FT=None))
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:DP:FT .:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=False)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=., DP=None, FT=None))
    >>> # missing call when missing calls excluded
    >>> _test_call_generator("chrom 1 . A T . PASS . GT:DP:FT ./.:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=True)

    >>> _test_call_generator("chrom 1 . A T . PASS . GT:DP:FT .:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=True)

    >>> # mix of snps and indels in same record, all calls are homozygous
    >>> _test_call_generator("chrom 1 . T G,TGA . PASS . GT:FT 1/1:PASS 2/2:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE1, CallData(GT=1/1, FT=PASS))
    >>> _test_call_generator("chrom 1 . T G,TGA . PASS . GT:FT 1/1:PASS 2/2:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE2, CallData(GT=2/2, FT=PASS))
    >>> _test_call_generator("chrom 1 . T G,TGA . PASS . GT:FT 1/1:PASS 2/2:PASS", exclude_snps=False, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE1, CallData(GT=1/1, FT=PASS))
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE2, CallData(GT=2/2, FT=PASS))
    >>> _test_call_generator("chrom 1 . T G,TGA . PASS . GT:FT 2/2:PASS 1/1:PASS", exclude_snps=False, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE1, CallData(GT=2/2, FT=PASS))
    Record(CHROM=chrom, POS=1, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE2, CallData(GT=1/1, FT=PASS))

    >>> # mix of heterozygous snps and indels in same call
    >>> _test_call_generator("chrom 2 . T G,TGA . PASS . GT:FT 1/2:PASS", exclude_snps=False, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=2, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE, CallData(GT=1/2, FT=PASS))
    >>> _test_call_generator("chrom 2 . T G,TGA . PASS . GT:FT 1/2:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=2, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE, CallData(GT=1/2, FT=PASS))
    >>> _test_call_generator("chrom 2 . T G,TGA . PASS . GT:FT 1/2:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=2, REF=T, ALT=[G, TGA]) Call(sample=SAMPLE, CallData(GT=1/2, FT=PASS))
    >>> _test_call_generator("chrom 2 . T G,TGA . PASS . GT:FT 1/2:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)

    """
    reader = vcf.VCFReader(input)
    for record in reader:
        for call in record.samples:
            # All calls are included by default, unless explicitly excluded
            # A call with both snp and indel is not excluded unless both snps and indels are excluded
            keep_snp_and_indel_mix = not (exclude_snps and exclude_indels) and is_snp_call(record, call) and is_indel_call(record, call)
            if not keep_snp_and_indel_mix and exclude_snps and is_snp_call(record, call):
                continue
            if not keep_snp_and_indel_mix and exclude_indels and is_indel_call(record, call):
                continue
            if exclude_vars and is_other_variant_call(record, call):
                continue
            if exclude_refs and is_ref_call(record, call):
                continue
            if exclude_hetero and call.is_het:
                continue
            if exclude_filtered and is_filtered_call(record, call):
                continue
            if exclude_missing and is_missing_call(record, call):
                continue

            yield (record, call)
