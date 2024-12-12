#check if other species behave the same in this sex-specific loci

## after generate call variants and generate vcf file of each species.
## extract targeted loci
vcftools --vcf /crex/proj/naiss2023-6-57/4-vcf/3-filtered_vcf_diploids/s.viminalis.snps.filter.vcf.gz --positions sher.allele_diff_sex.positions --recode --out positions-svim.snps.filter
gatk VariantsToTable -F CHROM -F POS -F TYPE -GF GT -GF AD -V positions-svim.snps.filter.recode.vcf -O positions-svim.snps.filter.gt.tbl

vcftools --gzvcf /crex/proj/naiss2023-6-57/4-vcf/3-filtered_vcf_diploids/s.viminalis.snps.filter.vcf.gz --positions sher.fheter-mhomo.positions --recode --out fheter-mhomo-svim.snps.filter
gatk VariantsToTable -F CHROM -F POS -F TYPE -GF GT -GF AD -V fheter-mhomo-svim.snps.filter.recode.vcf -O fheter-mhomo-svim.snps.filter.gt.tbl


for i in *pos.gt.tbl;do awk '{print $1,$2,$3,$4,$5,$10,$11,$12,$13,$6,$7,$8,$9,$14,$15}' $i | sed 's/ /\t/g' > ${i}"2"; done

perl allele_fdp_m0.pl spur-ssuc-6ind.indel.h1mask_She1.allele_fdp_m0.fhomo.pos.gt.tbl2 > spur-ssuc-6ind.indel.h1mask_She1.allele_fdp_m0.fhomo.pos.gt.tbl2.check
perl allele_fheter_mhomo.pl spur-ssuc-6ind.indel.h1mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2 > spur-ssuc-6ind.indel.h1mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2.check
perl allele_fheter_mhomo.pl  spur-ssuc-6ind.indel.h2mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2 > spur-ssuc-6ind.indel.h2mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2.check
perl allele_fdp_m0.pl  spur-ssuc-6ind.snps.h1mask_She1.allele_fdp_m0.fhomo.pos.gt.tbl2 > spur-ssuc-6ind.snps.h1mask_She1.allele_fdp_m0.fhomo.pos.gt.tbl2.check
perl allele_fheter_mhomo.pl  spur-ssuc-6ind.snps.h1mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2 > spur-ssuc-6ind.snps.h1mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2.check
perl allele_fheter_mhomo.pl  spur-ssuc-6ind.snps.h2mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2 > spur-ssuc-6ind.snps.h2mask_She1.allele_fheter_mhomo.strict.pos.gt.tbl2.check

