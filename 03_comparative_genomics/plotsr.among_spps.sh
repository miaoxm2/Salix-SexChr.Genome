### create synteny alignment between genome assemblies
### detect structural variation

#1. aligne the genome-minimap"
pwd="/crex/proj/snic2021-6-33/0xiaomeng/12-synteny/1-plotsr/pre-assembly/"

module load bioinfo-tools minimap2/2.24-r1122
module load bioinfo-tools samtools/1.17

minimap2 -ax asm5 -t 10 --eqx /crex/proj/snic2021-6-33/0xiaomeng/12-synteny/1-plotsr/0-pre-assembly/sarb.chr15.fa /crex/proj/snic2021-6-33/0xiaomeng/12-synteny/1-plotsr/0-pre-assembly/spur.chr15z.fasta  | samtools sort -@ 10 -O BAM - > sarb_spur-z.chr15.bam 
samtools index sarb_spur-z.chr15.bam 

minimap2 -ax asm5 -t 10 --eqx /crex/proj/snic2021-6-33/0xiaomeng/12-synteny/1-plotsr/0-pre-assembly/spur.chr15w.fasta /crex/proj/snic2021-6-33/0xiaomeng/12-synteny/1-plotsr/0-pre-assembly/svim.chr15w.fa | samtools sort -@ 10 -O BAM - > spur_svim.chr15w.bam 
samtools index spur_svim.chr15w.bam 

#step2. structural annotation"
module load bioinfo-tools MUMmer/3.23

syri -c 1-bam/sarb_spur-z.chr15.bam -r ../0-pre-assembly/sarb.chr15.fa -q ../0-pre-assembly/spur.chr15z.fasta --dir 2-syri/ --prefix sarb_spur-z.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

syri -c 1-bam/ssuc_spur.chr15z.bam -r ../0-pre-assembly/ssuc.chr15z.fa -q ../0-pre-assembly/spur.chr15z.fasta --dir 2-syri/ --prefix ssuc_spur.chr15.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

syri -c 1-bam/spur_sher.chr15z.bam -r ../0-pre-assembly/spur.chr15z.fasta -q ../0-pre-assembly/sher.h2.chr15.fasta --dir 2-syri/ --prefix spur_sher.chr15z.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

syri -c 1-bam/h2_h1.chr15.bam -r ../0-pre-assembly/sher.h2.chr15.fasta -q ../0-pre-assembly/sher.h1.chr15.fasta --dir 2-syri/ --prefix h2-h1.chr15.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

#syri -c 1-bam/h1_spur.chr15w.bam -r ../0-pre-assembly/sher.h1.chr15.fasta -q ../0-pre-assembly/spur.chr15w.fasta --dir 2-syri/ --prefix h1_spur.chr15w.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

syri -c 1-bam/spur_svim.chr15w.bam -r ../0-pre-assembly/spur.chr15w.fasta -q ../0-pre-assembly/svim.chr15w.fa --dir 2-syri/ --prefix spur_svim.chr15w.strict. -F B --invgaplen 1000000  --tdmaxolp 0.5 --hdrseq --cigar --maxsize 100000

#step3. check inversion. criteria: coverage>0.1
for file in 2-syri/*syri.out  
do
echo $file >> 2-syri-inv-check.out
awk '$11=="INVAL"' $file |awk '{a[$10]+=$3-$2;b[$10] += $7-$8;if (!min_ref_start[$10] || $2 < min_ref_start[$10]) min_ref_start[$10] = $2;if (!max_ref_end[$10] || $3 > max_ref_end[$10]) max_ref_end[$10] = $3;if (!max_query_start[$10] || $4 > max_query_start[$10]) max_query_start[$10] = $7;if (!min_query_end[$10] || $5 < min_query_end[$10]) min_query_end[$10] = $8;} END{for( i in a) print i,a[i],min_ref_start[i],max_ref_end[i],a[i]/(max_ref_end[i]-min_ref_start[i]),b[i],min_query_end[i],max_query_start[i],b[i]/(max_query_start[i]-min_query_end[i])}' >> 2-syri-inv-check.out
awk '$11=="INVAL"' $file | awk '{a[$10]+=$3-$2;b[$10] += $7-$8;if (!min_ref_start[$10] || $2 < min_ref_start[$10]) min_ref_start[$10] = $2;if (!max_ref_end[$10] || $3 > max_ref_end[$10]) max_ref_end[$10] = $3;if (!max_query_start[$10] || $4 > max_query_start[$10]) max_query_start[$10] = $7;if (!min_query_end[$10] || $5 < min_query_end[$10]) min_query_end[$10] = $8;} END{for( i in a) if(a[i]/(max_ref_end[i]-min_ref_start[i])<0.1 || b[i]/(max_query_start[i]-min_query_end[i])<0.1) print i}' >> 2-syri-inv-check.out
done

#step4. filter inversion according to the output of last step
#h2-h2 delete INV27 INV29 INV25
grep 'INV27' 2-syri/h2-h1.chr15.strict.syri.out |awk '$11=="INVAL"' | awk '{print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,"INV",$12}' | sed 's/ /\t/g' > 2-syri/h2-h1.chr15.strict.tmpinv27
grep 'INV29' 2-syri/h2-h1.chr15.strict.syri.out |awk '$11=="INVAL"' | awk '{print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,"INV",$12}' | sed 's/ /\t/g' > 2-syri/h2-h1.chr15.strict.tmpinv29
grep 'INV25' 2-syri/h2-h1.chr15.strict.syri.out |awk '$11=="INVAL"' | awk '{print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,"INV",$12}' | sed 's/ /\t/g' > 2-syri/h2-h1.chr15.strict.tmpinv25
grep -v 'INV27' 2-syri/h2-h1.chr15.strict.syri.out |grep -v 'INV25' | grep -v 'INV29' > 2-syri/h2-h1.chr15.strict.tmpall
cat 2-syri/h2-h1.chr15.strict.tmp* > 2-syri/h2-h1.chr15.strict.syri.filter.out

#sarb-spur-z delete INV27 INV23
grep 'INV27' 2-syri/sarb_spur-z.strict.syri.out |awk '$11=="INVAL"' | awk '{print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,"INV",$12}' | sed 's/ /\t/g' > 2-syri/sarb_spur-z.strict.tmpinv27
grep 'INV23' 2-syri/sarb_spur-z.strict.syri.out |awk '$11=="INVAL"' | awk '{print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,"INV",$12}' | sed 's/ /\t/g' > 2-syri/sarb_spur-z.strict.tmpinv23
grep -v 'INV27' 2-syri/sarb_spur-z.strict.syri.out |grep -v 'INV23' > 2-syri/sarb_spur-z.strict.tmpall
cat 2-syri/sarb_spur-z.strict.tmp* > 2-syri/sarb_spur-z.strict.syri.filter.out

#step5. plot
plotsr --sr 2-syri/sarb_spur-z.strict.syri.filter.out --sr 2-syri/spur_sher.chr15z.strict.syri.out --sr 2-syri/h2-h1.chr15.strict.syri.filter.out --sr 2-syri/h1_spur.chr15w.strict.syri.out --sr 2-syri/spur_svim.chr15w.strict.syri.out --genomes genomes.sarb-spur-sher-svim.txt -o plotsr.chr15.sarb-spur-sher-svim.pdf --cfg base.cfg -H 6 -W 6 -S 0.7 -R --markers markers.bed

