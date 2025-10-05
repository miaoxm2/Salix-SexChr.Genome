## calculate fst between sexes
module load bioinfo-tools vcftools/0.1.16
#fst she1 f-m
vcftools --vcf h1mask_She1-20ind.snps.filter.vcf \
--weir-fst-pop 4-fst/she1.f10.pop --weir-fst-pop 4-fst/she1.m10.pop \
--max-meanDP 200 --min-meanDP 10  --maf 0.1 --fst-window-size 10000 --fst-window-step 1000 \
--out 4-fst/h1mask_She1-20ind.f-m.10k2

vcftools --vcf h2mask_She1-20ind.snps.filter.vcf \
--weir-fst-pop 4-fst/she1.f10.pop --weir-fst-pop 4-fst/she1.m10.pop \
--max-meanDP 200 --min-meanDP 10  --maf 0.1 --fst-window-size 10000 --fst-window-step 1000 \
--out 4-fst/h2mask_She1-20ind.f-m.10k
