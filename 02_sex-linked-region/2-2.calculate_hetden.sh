## count het. number per ind. 
# split SNP file into several files. Each chr per file
awk '/^h1_chr/{print > $1".She1-20ind.snp.vcf"}' ../h1mask_She1-20ind.snps.filter.vcf

#count the number of hete site in window
for file in *vcf; do perl ind-heter_rate.per_chr.10kb.pl $file > $file.hetden10kb; done

for file in *vcf.hetden10kb; do
awk -F '\t' '
  function basename(file) {
    sub(".She1.*", "", file)
    return file
  }
  {female=($2+$3+$4+$5+$6+$7+$8+$9+$10+$11)/10;male=($12+$13+$14+$15+$16+$17+$18+$19+$20+$21)/10;print basename(FILENAME),$1,female,male}' $file > $file.ave
done

for file in *ave; do tail -n +2 $file >> h1_mask.She1-20ind.snp.vcf.hetden10kb.ave; done
