### filter variants

module load bioinfo-tools GATK/4.3.0.0

ls $indir/$prefix*g.vcf.gz > $indir/$prefix.g.vcf.list

gatk CombineGVCFs -R $ref -V $indir/$prefix.g.vcf.list -O $outdir/$prefix.g.vcf.gz

gatk GenotypeGVCFs -R $ref -V $outdir/$prefix.g.vcf.gz -O $outdir/$prefix.raw.vcf.gz

gatk SelectVariants -R $ref -V $outdir/$prefix.raw.vcf.gz -O $outdir/$prefix.snps.raw.vcf -select-type SNP

gatk SelectVariants -R $ref -V $outdir/$prefix.raw.vcf.gz -O $outdir/$prefix.indel.raw.vcf -select-type INDEL

## 1-hard filtering pipelines

#indels
gatk  --java-options "-Xmx30G" VariantFiltration \
    -V $prefix.indel.raw.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -filter "DP < 30" --filter-name "DP30" \
    -O $prefix.indel.filter1.vcf

#snp
gatk  --java-options "-Xmx30G" VariantFiltration \
    -V $prefix.snps.raw.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "DP < 30" --filter-name "DP30" \
    -O $prefix.snps.filter1.vcf

## 2-filter for downstream analysis

#filter vcf -biallelic and remove marked
#vcftools --vcf $vcf --out $filter --recode --remove-filtered-all --min-alleles 2 --max-alleles 2 

#filter vcf -depth and quality 
MAF=0.1
MISS=1 # --max-missing $MISS 
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=100

vcftools --vcf $filter.recode.vcf --out $filter2 --recode --remove-filtered-all --maf $MAF --minQ $QUAL --minDP $MIN_DEPTH --maxDP $MAX_DEPTH


