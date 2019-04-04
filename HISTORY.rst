.. :changelog:

History
=======

1.2.0 (2019-04-04)
---------------------
* Fix defect in narrow command wrongly printing ALT=. when GT=.
* Add the ``count`` command to count samples, positions, calls, snps, indels,
  other variants, filtered calls, missing calls, and filter reasons.
* Add the ``plot`` command to plot calls along the length of the genome and show
  the location of filtered calls.
* Change the text of the compare report to refer to "Calls", not "Sample snps".
* Drop support for Python 3.4, which is not supported by matplotlib.
* Add support for Python 3.7.

1.1.1 (2019-03-26)
---------------------
* Replace None with '.' when printing call data.
* Support VCF files with multiple alternate alleles per position.

1.1.0 (2019-02-06)
---------------------
* Support reading gzip compressed vcf files.


1.0.0 (2018-11-20)
---------------------

* First public release.
