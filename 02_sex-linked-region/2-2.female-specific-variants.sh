##find out different alleles among all female/males

perl allele_diff_sex.pl $filter2.bi.gt.tbl > $filter2.bi.gt.allele_diff_sex
perl allele_f0_mdp.pl $filter2.bi.gt.tbl > $filter2.bi.gt.allele_f0_mdp
perl allele_fdp_m0.pl $filter2.bi.gt.tbl > $filter2.bi.gt.allele_fdp_m0
perl allele_fheter_mhomo.pl $filter2.bi.gt.tbl > $filter2.bi.gt.allele_fheter_mhomo
perl allele_fhomo_mheter.pl $filter2.bi.gt.tbl > $filter2.bi.gt.allele_fhomo_mheter


awk '{print $1,"\t",$2}' $filter2.bi.gt.allele_diff_sex >$filter2.bi.gt.allele_diff_sex.positions

vcftools --vcf $filter2.bi.recode.vcf \
--positions $filter2.bi.gt.allele_diff_sex.positions \
--recode --out $filter2.bi.gt.allele_diff_sex

# extract the gene annotation for these alleles
gff="/crex/proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/Sherbacea_final.h1.gene.gff"
module load bioinfo-tools BEDTools/2.31.1
bedtools intersect -a $filter2.bi.gt.allele_diff_sex.recode.vcf -b $gff -wa -wb > $filter2.bi.gt.allele_diff_sex.gene.vcf
