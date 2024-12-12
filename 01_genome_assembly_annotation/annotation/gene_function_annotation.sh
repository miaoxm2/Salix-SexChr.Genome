### gene function annotation

## interproscan
module load bioinfo-tools InterProScan/5.62-94.0
interproscan.sh -i h1_final.braker_combined.try1.renamed.fixed_cleaned.pro.fa -f GFF3 -o h1_final.interproscan_output.gff --goterms --pathways --iprlookup

## Eggnog-mapper. use online version
http://eggnog-mapper.embl.de/

## TAIR11
module load bioinfo-tools blast
#makeblastdb -in araport11.pro.fasta -dbtype prot -title database_araport11 -parse_seqids -out araport11.pro.db
blastp -query ../h1_final.braker_combined.try1.renamed.fixed_cleaned.pro.fa -out araport11-h1.out -db araport11.pro.db -outfmt 6 -evalue 1e-5 -num_threads 8 -max_target_seqs 1

blastp -query ../../h2/5-combine/h2_final.braker_combined.try1.renamed.fixed.pro.fa  -out araport11-h2.out -db araport11.pro.db -outfmt 6 -evalue 1e-5 -num_threads 8 -max_target_seqs 1

