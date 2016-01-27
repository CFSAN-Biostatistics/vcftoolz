What is in this directory
-------------------------
compare-vcf.sh : uses bcftools isec to compare VCF files
comparevcf.py  : custom python code to compare VCF files
preprocess.sh  : script to fix issues in the VCF files that prevent comparisons

Issues with VCF Files
---------------------
1. The chrom names are different in lyveset and snp-pipeline VCF files
   Resolved by sed in preprocess.sh

2. PyVCF thinks the lyveset ##contig headers are malformed
   Resolved by sed in preprocess.sh

3. PyVCF thinks some of the lyveset VCF records are not snp records.
   Added -a option to process all VCF records.
   Q: How does PyVCF determine if a record is a snp record?
   A: All the ALT alleles must be in this set: ['A', 'C', 'G', 'T']. 
      However, the newer version of PyVCF on github requires all the alleles
      to be in this set: ['A', 'C', 'G', 'T', 'N', '*']

4. Lyveset sometimes encodes alt allele=N.
   TODO: Determine whether this is a failed snp.
   Added check to ignore base==N when collecting the list of snps.

5. Snp-pipeline might not be encoding failed snps genotype (GT) properly.
   Snp-pipeline possible GT values:
      GT=0 : most common base matches the reference
      GT=1 : snp with allele matching the first allele ALT column
      GT=. : no good depth. unable to make a call
   Snp-pipeline possible ALT values:
      ALT=.      : either no good depth or the most common base matches the reference
      ALT=[ACTG] : most common base when it does not match the reference


VCF Spec
--------
If a call cannot be made for a sample at a given locus, "."  should be specified for each missing allele in the GT field.
