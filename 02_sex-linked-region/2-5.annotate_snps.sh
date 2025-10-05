## filter vcf files
module load bioinfo-tools GATK/4.3.0.0
module load bioinfo-tools vcftools/0.1.16
##snps- after hard filter
filter1="h1mask_She1-20ind.snps.filter.vcf"
filter2="h1mask_She1-20ind.snps.filter2"

gatk VariantsToTable -F CHROM -F POS -F TYPE -GF GT -V $filter1 -O $filter1.gt.tbl

#filter for deeper
### set filters
MAF=0.1
MISS=0.8
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=100

vcftools --vcf $filter1 \
--maf $MAF  --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--recode --remove-filtered-all \
--min-alleles 2 --max-alleles 2 \
--out $filter2.bi

## annotate snps
module load bioinfo-tools snpEff/5.2
vcf="h1mask_She1-20ind.snps.filter2.bi.recode.vcf "

snpEff -c snpEff.config -formatEff  -stats h1mask_She1-20ind.eff -v sher.h1 $vcf > h1mask_She1-20ind.snp.eff.vcf
grep -v '#' h1mask_She1-20ind.snp.eff.vcf | awk '$7=="PASS"''{split($8, a, "|");print $1,$2,a[1],a[2],a[6]}' > h1mask_She1-20ind.snp.eff.function
