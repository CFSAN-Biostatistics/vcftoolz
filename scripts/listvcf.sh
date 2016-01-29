# List a VCF file or just one position from a VCF file in nicely formatted columns.
#
# Examples
# --------
# ./listvcf.sh file.vcf
#
# ./listvcf.sh file.vcf  2345 # just position 2345


if [[ -n $2 ]]; then
    grep CHROM $1 > $1.qry
    grep $2 $1 >> $1.qry
    column -t $1.qry | less -S
else
    cat $1 | grep -v "##" | column -t  | less -N -S
fi
