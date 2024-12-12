## detect hemizygous sites

#1- calculate depth

for bam in h2mask_She2-*markdup.addrg.bam
do
        prefix="coverage/"${bam%%\.*}
        mosdepth --by 10000 --mapq 20 --threads 5 --fast-mode $out $bam
done

cd coverage/
for file in h2mask_She2-*regions.bed.gz
do
echo ${file%%.*} >${file%%.*}.ave 
awk 'NR==FNR{ total += $4;len=NR } NR>FNR{ print  $4 /(total/len)}' <(zcat $file ) <(zcat $file) >> ${file%%.*}.ave 
done

depth="/crex/proj/snic2021-6-33/vinciane/dwarf/2-map_bam/coverage/h1mask_She.32ind.10kb.ave-depth"
awk -v OFS='\t' '{for (i=23;i<=28;i++){f=f+$i};for (i=29;i<=34;i++){m=m+$i};if (f+m>0){ave=f/(f+m)};print $1,$2,f,m,ave; f=0;m=0}' $depth > h1mask_She2.f-m.10k.depth


#2- generate dp<5 bed file
mkdir per-ind
for input in h1mask_She2*per-base.bed.gz
do
name=${input##*/}
sp=${name%.per-base.bed.gz}
zcat $input | awk '$4<5' > per-ind/$sp.dp0
done

for file in per-ind/h1mask_She2-*.dp0; do awk '($1=="h1_chr15"){a=a+$3-$2;}END{print a}' $file; echo $file; echo; done


#3- intersect bed file
module load bioinfo-tools BEDTools/2.31.1

#find all F/M with dp<5
multiIntersectBed -i  h1mask_She2-M*.dp0 > h1mask_She2-M.dp0

awk '$4==10' h1mask_She2-M.dp0 > h1mask_She2-M.dp0.all

multiIntersectBed -i  per-ind/h1mask_She2-F*.dp0 > h1mask_She2-F.dp0

awk '$4==10' h1mask_She2-F.dp0 > h1mask_She2-F.dp0.all

## Conditions:
# M=0 & F!=0
bedtools intersect -a h1mask_She2-M.dp0.all -b h1mask_She2-F.dp0.all -v > h1mask_She2-dp.M0F1
cut -f1,2,3 h1mask_She2-dp.M0F1 > h1mask_She2-dp.M0F1.bed

# F=0 & M!=0
bedtools intersect -a h1mask_She2-F.dp0.all -b h1mask_She2-M.dp0.all -v > h1mask_She2-dp.F0M1
cut -f1,2,3 h1mask_She2-dp.F0M1 > h1mask_She2-dp.F0M1.bed


# F=0 & M=0
bedtools intersect -a h1mask_She2-F.dp0.all -b h1mask_She2-M.dp0.all -wa -wb > h1mask_She2-dp.F0M0
cut -f1,2,3 h1mask_She2-dp.F0M0 > h1mask_She2-dp.F0M0.bed

