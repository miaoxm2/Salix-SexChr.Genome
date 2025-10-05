#calculate depth per individual per sites.
#same pipeline applied to h2

##generate depth using mosdepth
for bam in h1mask_She1-*markdup.addrg.bam
do
        prefix="coverage/"${bam%%\.*}
        mosdepth --by 10000 --mapq 20 --threads 5 --fast-mode $prefix $bam
done

##merge depth into list 
cut -f1-2 <(zcat /crex/proj/naiss2023-6-57/3-markdup/coverage/h1mask_She1-F05.ave) > /crex/proj/naiss2023-6-57/3-markdup/coverage/h1mask_She1.coord
paste h1mask_She1.coord /crex/proj/naiss2023-6-57/3-markdup/coverage/h1mask_She1*ave > h1mask_She1.10kb.ave-depth

##find out sites where depth<5
for input in /crex/proj/naiss2023-6-57/3-markdup/coverage/*per-base.bed.gz
do
name=${input##*/}
sp="per-ind/"${name%.per-base.bed.gz}
zcat $input | awk '$4<=5' > $sp.dp0
done

##calculate average depth
cd coverage/
for file in h1mask_She1-*regions.bed.gz
do
echo ${file%%.*} >${file%%.*}.ave
awk 'NR==FNR{ total += $4;len=NR } NR>FNR{ print  $4 /(total/len)}' <(zcat $file ) <(zcat $file) >> ${file%%.*}.ave
done

cut -f1-2 <(zcat h1mask_She1-F05.regions.bed.gz) > coord
paste coord *ave > h1mask_She.20ind.10kb.ave-depth


#create normalised depth from mosdepth output 
for file in h1mask_She1*regions.bed.gz; do awk 'NR==FNR{ total += $4;len=NR } NR>FNR{ print  $1,$2,$4 /(total/len)}' <(zcat $file ) <(zcat $file) > ${file%%.*}.ave ; done
##create list
ls h1mask_She1*ave | sed 's/\t//g' > h1mask_She1.list

##merge *ave together
perl join_dp.pl h1mask_She1.list > h1mask_She1.20ind.ave
