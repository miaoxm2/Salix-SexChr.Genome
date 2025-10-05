#### find homologs between z and w chromosome

module load bioinfo-tools blast/2.14.1+
#1-soft blast -- to shrink the range of specific genes
#evalue<1e5; identity>60; max target<5; alignment/query length>50%
blastp -query ../library/sher.h1_chr15.pro.fasta -db ../library/sher.h2_chr15.pro.db -out sher.h1chr15-h2chr15.pro.soft.blast -evalue 1e-5 -max_target_seqs 5 -num_threads 5 -outfmt 6 -qcov_hsp_perc 50
blastp -query ../library/sher.h2_chr15.pro.fasta -db ../library/sher.h1_chr15.pro.db -out sher.h2chr15-h1chr15.pro.soft.blast -evalue 1e-5 -max_target_seqs 5 -num_threads 5 -outfmt 6 -qcov_hsp_perc 50
awk '$3>60' sher.h1chr15-h2chr15.pro.soft.blast | awk '{print $1}' |sed 's/ //g' >sher.h1chr15-h2chr15.pro.soft.id60.h1genes
awk '$3>60' sher.h1chr15-h2chr15.pro.soft.blast | awk '{print $2}' |sed 's/ //g' >sher.h1chr15-h2chr15.pro.soft.id60.h2genes
awk '$3>60' sher.h2chr15-h1chr15.pro.soft.blast | awk '{print $2}' |sed 's/ //g' >sher.h2chr15-h1chr15.pro.soft.id60.h1genes
awk '$3>60' sher.h2chr15-h1chr15.pro.soft.blast | awk '{print $1}' |sed 's/ //g' >sher.h2chr15-h1chr15.pro.soft.id60.h2genes


#2-stringent blast -- to focus on shared genes
#evalue<1e10; identity>80; max target<5; alignment/query length>70%
blastp -query ../library/sher.h1_chr15.pro.fasta -db ../library/sher.h2_chr15.pro.db -out sher.h1chr15-h2chr15.pro.hard.blast -evalue 1e-10 -max_target_seqs 5 -num_threads 5 -outfmt 6 -qcov_hsp_perc 70

blastp -query ../library/sher.h2_chr15.pro.fasta -db ../library/sher.h1_chr15.pro.db -out sher.h2chr15-h1chr15.pro.hard.blast -evalue 1e-10 -max_target_seqs 5 -num_threads 5 -outfmt 6 -qcov_hsp_perc 70

awk '$3>80' sher.h1chr15-h2chr15.pro.hard.blast | awk '{print $1,"=",$2}' |sed 's/ //g' >sher.h1chr15-h2chr15.pro.hard.id80.pairs
awk '$3>80' sher.h2chr15-h1chr15.pro.out | awk '{print $2,"=",$1}' |sed 's/ //g' >sher.h2chr15-h1chr15.pro.hard.id80.pairs

cat sher.h1chr15-h2chr15.pro.soft.id60.h1genes sher.h2chr15-h1chr15.pro.soft.id60.h1genes | sort | uniq > sher.soft.shared.h1genes
cat sher.h1chr15-h2chr15.pro.soft.id60.h2genes sher.h2chr15-h1chr15.pro.soft.id60.h2genes | sort | uniq > sher.soft.shared.h2genes

awk 'FNR==NR { gene_ids[$1]++; next } !($1 in gene_ids)' sher.soft.shared.h1genes ../library/sher.h1_chr15.gene.list >sher.soft.only.h1genes
awk 'FNR==NR { gene_ids[$1]++; next } !($1 in gene_ids)' sher.soft.shared.h2genes ../library/sher.h2_chr15.gene.list >sher.soft.only.h2genes


#-outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore ppos qcovs qcovhsp"
#-max_target_seqs <Integer, >=1> Maximum number of aligned sequences to keep 
#-evalue <Real>   Expectation value (E) threshold for saving hits.
#-qcov_hsp_perc <Real, 0..100>   Percent query coverage per hsp


### detect pseudo genes 

## 1- using repeat&cds-masked genome as query database
qry="sher.soft.only.h1pro.fasta"
db="/crex/proj/snic2021-6-33/0xiaomeng/14-geneloss/library/genome.chr15/sher.chr15z.mask_rep_cds.db"
pre="sher.onlyh1-h2mask"
tblastn -query $qry -db $db -out $pre.tblastn -evalue 1e-5 -num_threads 5 -outfmt 6 

qry="sher.soft.only.h2pro.fasta"
db="/crex/proj/snic2021-6-33/0xiaomeng/14-geneloss/library/genome.chr15/sher.chr15w.mask_rep_cds.db"
pre="sher.onlyh2-h1mask"
tblastn -query $qry -db $db -out $pre.tblastn -evalue 1e-5 -num_threads 5 -outfmt 6 

qry="sher.soft.only.h1pro.fasta"
db="/crex/proj/snic2021-6-33/0xiaomeng/14-geneloss/library/genome.chr15/sher.chr15w.mask_rep_cds.db"
pre="sher.onlyh1-h1mask"
tblastn -query $qry -db $db -out $pre.tblastn -evalue 1e-5 -num_threads 5 -outfmt 6 

qry="sher.soft.only.h2pro.fasta"
db="/crex/proj/snic2021-6-33/0xiaomeng/14-geneloss/library/genome.chr15/sher.chr15z.mask_rep_cds.db"
pre="sher.onlyh2-h2mask"
tblastn -query $qry -db $db -out $pre.tblastn -evalue 1e-5 -num_threads 5 -outfmt 6 


## 2- generate bed file for pseudo and putative gene loss. use h1-h1mask as an example
awk '$3>40 && $4>30 {print $1,$9,$10}' sher.h1-h1mask.tblastn  > sher.h1-h1mask.tblastn.id40len30.bed

perl 1-sort_gene.merge.pl sher.h1-h1mask.tblastn.id40len30.bed > sher.h1-h1mask.tblastn.id40len30.gene
perl 2-sort_coord.merge.pl sher.onlyh1-h1mask.id40len30.bed | awk '{print "h1_chr15",$0}' | sed 's/ /\t/g' > sher.onlyh1-h1mask.pseudo.bed 

