===============================
VCF Toolz
===============================


.. Image showing the PyPI version badge - links to PyPI
.. image:: https://img.shields.io/pypi/v/vcftoolz.svg
        :target: https://pypi.python.org/pypi/vcftoolz

.. Image showing the Travis Continuous Integration test status, commented out for now
.. .. image:: https://img.shields.io/travis/CFSAN-Biostatistics/vcftoolz.svg
..        :target: https://travis-ci.org/CFSAN-Biostatistics/vcftoolz

.. Image showing the JOSS paper badge
.. image:: http://joss.theoj.org/papers/10.21105/joss.01144/status.svg
   :target: https://doi.org/10.21105/joss.01144

Tools for working with Variant Call Format files.

VCF Toolz was developed by the United States Food
and Drug Administration, Center for Food Safety and Applied Nutrition.

* Free software
* Documentation: https://vcftoolz.readthedocs.io
* Source Code: https://github.com/CFSAN-Biostatistics/vcftoolz
* PyPI Distribution: https://pypi.python.org/pypi/vcftoolz


Features
--------

* Compares the snps in two or more VCF files.
* Lists the snps that are unique to each VCF file with full genotype information per snp.
* Lists the snps that are missing from each VCF file if present in at least two other VCF files.
* Generates Venn diagrams of positions and snps in the VCF files.
* Reports precision, recall, and F1 score when the truth is known.
* Reports the effectiveness of filtered variants when the truth is known.
* Reformat the VCF file in a tall-narrow format for easy viewing and diffs.
* Count samples, positions, calls, snps, indels, other variants, missing calls, and filter reasons.
* Plot calls along the length of the genome and show the location of filtered calls.


Citing VCF Toolz
--------------------------------------

To cite VCF Toolz, please reference the VCF Toolz paper:

    https://doi.org/10.21105/joss.01144


License
-------

See the LICENSE file included in the VCF Toolz distribution.

