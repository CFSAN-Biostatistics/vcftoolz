---
title: 'vcftoolz: a Python package for comparing and evaluating Variant Call Format files.'
tags:
  - bioinformatics
  - next-generation sequencing
  - variant call format
  - VCF
authors:
 - name: Steve Davis
   orcid: 0000-0002-5581-1922
   affiliation: "1"
affiliations:
 - name: U.S. Food and Drug Administration
   index: 1
date: 21 November 2018
bibliography: paper.bib
---

# Summary and Need Statement

The analysis of next-generation sequence data often involves variant calling --
a process for discovering genomic variants, which are differences between genomes.
The standard output file format of variant callers is the Variant Call Format
(VCF) (https://github.com/samtools/hts-specs).

Researchers have a need to view and compare VCF files.  The behaviour of
alternative variant calling algorithms can be analysed by comparing the VCF files
produced by those algorithms.  The performance of a variant calling algorithm
can be evaluated by comparing against a known truth VCF dataset.

Here, we present ``vcftoolz``, software to facilitate comparing and evaluating
the variant calls in VCF files.  The core functionality of ``vcftoolz`` is the
capability to compare two or more VCF files, producing a report, Venn Diagrams,
and a spreadsheet identifying the concordance between the VCF files.  The artifacts 
produced by ``vcftools`` are not available from other tools.

# Related Research

The ``vcftoolz`` software is being used as part of an ongoing effort to compare
and evaluate the variant callers used by multiple government agencies involved
in the analysis of pathogenic organisms of interest to food safety. In this effort,
we use multiple variant callers to construct VCF files from closely-related pathogens.
The ``vcftoolz`` software identifies the concordance between the VCF files produced
by the alternative variant callers and facilitates algorithm improvements.
 
# Prior Related Work

The ``RTG Tools`` package [@RTGtools] has advanced capabilities to compare VCF files containing complex
variants, but does not support VCF files with multiple samples per file.

The ``BCFtools`` [@BCFtools] package has the capability to create intersections, unions and complements of VCF files,
as well as other useful tools for working with VCF files.

The ``VCFtools`` [@VCFtools] package has the capability to calculate differences between VCF files, among other functions.

# Links

Documentation:
https://vcftoolz.readthedocs.io/en/latest/readme.html

Source Code:
https://github.com/CFSAN-Biostatistics/vcftoolz

PyPI Distribution:
https://pypi.python.org/pypi/vcftoolz


# References
