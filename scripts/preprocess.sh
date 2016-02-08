# Shorten snp-pipeline chrom names to match lyveset
sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006053_AOZS01000000_AOZS01000001-AOZS01000042/CP006053/g' snp-pipeline-0.4.2-bareilly.vcf
sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006054_AOZS01000000_AOZS01000001-AOZS01000042/CP006054/g' snp-pipeline-0.4.2-bareilly.vcf

sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006053_AOZS01000000_AOZS01000001-AOZS01000042/CP006053/g' snp-pipeline-0.5.1-bareilly.faildot.vcf
sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006054_AOZS01000000_AOZS01000001-AOZS01000042/CP006054/g' snp-pipeline-0.5.1-bareilly.faildot.vcf

sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006053_AOZS01000000_AOZS01000001-AOZS01000042/CP006053/g' snp-pipeline-0.5.1-bareilly.fail0.vcf
sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006054_AOZS01000000_AOZS01000001-AOZS01000042/CP006054/g' snp-pipeline-0.5.1-bareilly.fail0.vcf

sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006053_AOZS01000000_AOZS01000001-AOZS01000042/CP006053/g' snp-pipeline-0.5.1-bareilly.fail1.vcf
sed -i 's/Salmonella_enterica_subsp._enterica_serovar_Bareilly_str._CP006054_AOZS01000000_AOZS01000001-AOZS01000042/CP006054/g' snp-pipeline-0.5.1-bareilly.fail1.vcf

# Remove malformed contig lines from lyveset vcf file
sed -i '/##contig/d' lyveset-1.1.3-bareilly.vcf
sed -i '/##contig/d' lyveset-1.1.4-bareilly.withNs.vcf
sed -i '/##contig/d' lyveset-1.1.4-bareilly.noNs.vcf
