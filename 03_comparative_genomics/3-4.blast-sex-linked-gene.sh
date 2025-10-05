module load bioinfo-tools blast
#candidate gene
qry="/crex/proj/snic2021-6-33/0xiaomeng/8-blast/sher.candidate.pro.fasta"

#search in all available Salix and Populus genomes.

for file in /8-blast/library.pro/library_s*/*pdb
do
name=${file##*/}
db=${file%\.pdb}
sp=${name%.db.pdb}
blastp -query $qry -out candidate-$sp.out -db $db -outfmt 6 -evalue 1e-5 -num_threads 8
done
