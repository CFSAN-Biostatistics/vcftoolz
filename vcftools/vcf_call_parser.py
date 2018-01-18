#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
Module to parse vcf files and generate calls with their containing records.
The records and calls are in PyVCF format.
"""

from __future__ import print_function

import sys
import vcf

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
    """
    if len(record.REF) > 1:
        return False

    if is_ref_call(record, call):
        return False

    alleles = call.gt_bases.split(call.gt_phase_char())
    for bases in alleles:
        if len(bases) > 1:
            return False

    return True


def is_indel_call(record, call):
    """Return True if a sample call is an insertion or deletion.

    An indel is a call where ...

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is an indel.
    """
    return record.is_indel


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
    """
    if is_snp_call(record, call):
        return False
    if is_indel_call(record, call):
        return False
    if call.is_variant: # at least one GT allele is different from the REF
        # is_variant can be:
        # None  : gt == .
        # True  : gt > 0
        # False : gt == 0
        return True
    return False


def is_ref_call(record, call):
    """Return True if a sample call is a reference call.

    An ref call is a call where ...

    Parameters
    ----------
    record : pyvcf record
        Record parsed by PyVCF.
    call : pyvcf call
        Genotype call parsed by PyVCF.

    Returns
    -------
    True if the call is a reference call.
    """
    return call.gt_type == 0 # homogeneous ref


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
    """
    try:
        filt = call.data.FT     # If there are filters per sample, use them
    except AttributeError:
        filt = record.FILTER    # Otherwise, use the filter for the whole record

    if not filt or filt == "PASS" : # filter is not set or set to PASS
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
    """
    return call.gt_type is None and not any(call.data[1:])


def test_call_generator(test_record, exclude_snps=None, exclude_indels=None, exclude_vars=None, exclude_refs=None, exclude_hetero=None, exclude_filtered=None, exclude_missing=None):
    """This is a helper function used by the doctests to test various types of calls.
    """
    if sys.version_info < (3,):
        from StringIO import StringIO
    else:
        from io import StringIO
    in_memory_file = StringIO()
    in_memory_file.name = "test.vcf"
    print("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE", file=in_memory_file)
    print(test_record, file=in_memory_file)
    in_memory_file.seek(0)
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
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # snp when snps excluded
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # indel when snps included
    >>> test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # ref call when snps included
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # insertion when indels included
    >>> test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[AT]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # insertion when indels excluded
    >>> test_call_generator("chrom 1 . A AT . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # deletion when indels included
    >>> test_call_generator("chrom 1 . AT A . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=False, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=AT, ALT=[A]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # deletion when indels excluded
    >>> test_call_generator("chrom 1 . AT A . PASS . GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    
    >>> # structural variant when "other variants" included
    >>> test_call_generator("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=False, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=AT, ALT=[A]) Call(sample=SAMPLE, CallData(GT=1/1, FT=PASS))
    >>> # structural variant when "other variants" excluded
    >>> test_call_generator("chrom 1 . AT A . PASS SVTYPE=DEL; GT:FT 1/1:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    
    >>> # ref call when ref calls included
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=False, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=0/0, FT=PASS))
    >>> # ref call when ref calls excluded
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 0/0:PASS", exclude_snps=True, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # heterozygous snp call when snps and heterozygous calls included
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 0/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=False, exclude_filtered=True, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=0/1, FT=PASS))
    >>> # heterozygous snp call when snps included, but heterozygous calls excluded
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 0/1:PASS", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # FT filtered snp when snps and filtered calls included
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:FAIL", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=False, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1, FT=FAIL))
    >>> # FT filtered snp when snps included, but filtered calls excluded
    >>> test_call_generator("chrom 1 . A T . PASS . GT:FT 1/1:FAIL", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # FILTER filtered snp when snps and filtered calls included
    >>> test_call_generator("chrom 1 . A T . FAIL . GT 1/1", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=False, exclude_missing=True)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=1/1))
    >>> # FILTER filtered snp when snps included, but filtered calls excluded
    >>> test_call_generator("chrom 1 . A T . FAIL . GT 1/1", exclude_snps=False, exclude_indels=True, exclude_vars=True, exclude_refs=True, exclude_hetero=True, exclude_filtered=True, exclude_missing=True)

    >>> # missing call when missing calls included
    >>> test_call_generator("chrom 1 . A T . PASS . GT:DP:FT ./.:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=False)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=./., DP=None, FT=None))
    >>> test_call_generator("chrom 1 . A T . PASS . GT:DP:FT .:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=False)
    Record(CHROM=chrom, POS=1, REF=A, ALT=[T]) Call(sample=SAMPLE, CallData(GT=., DP=None, FT=None))
    >>> # missing call when missing calls excluded
    >>> test_call_generator("chrom 1 . A T . PASS . GT:DP:FT ./.:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=True)

    >>> test_call_generator("chrom 1 . A T . PASS . GT:DP:FT .:.:.", exclude_snps=False, exclude_indels=False, exclude_vars=False, exclude_refs=False, exclude_hetero=False, exclude_filtered=False, exclude_missing=True)

    """
    reader = vcf.VCFReader(input)
    for record in reader:
        for call in record.samples:
            # All calls are included by default, unless explicitly excluded
            #print(is_snp_call(record, call), call, file=sys.stdout)
            if exclude_snps and is_snp_call(record, call):
                continue
            if exclude_indels and is_indel_call(record, call):
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
