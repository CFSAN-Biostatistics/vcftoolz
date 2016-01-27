# This script runs bcftools isec to compare two vcf files.
# it compresses and indexes the vcf files as needed for isec.
#
# Example:
#    compare-vcf.sh  file1.vcf  file2.vcf

file1=$1
file2=$2

echo Zip vcf files
bgzip -c $file1 > $file1.gz
bgzip -c $file2 > $file2.gz

echo Index the zipped vcf files
tabix -f -p vcf $file1.gz
tabix -f -p vcf $file2.gz

echo Find intersections and differences
bcftools isec  --collapse none  -n=1 -p isec-out  $file1.gz  $file2.gz
