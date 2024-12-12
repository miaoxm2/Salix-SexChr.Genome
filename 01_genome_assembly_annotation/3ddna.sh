#! /bin/bash -l 
#SBATCH -A naiss2023-5-103  
#SBATCH -p core 
#SBATCH -n 5
#SBATCH -t 2-00:00:00 
#SBATCH -J 3ddna
#SBATCH --mail-user xiaomeng.mao@ebc.uu.se
#SBATCH --mail-type=ALL
module load bioinfo-tools lastz/1.04.00
#rm -f Sherbacea_v2.h1.*.cprops Sherbacea_v2.h1.mnd.*.txt
run-asm-pipeline.sh -i 50000 ../juicer/reference/Sherbacea_v2.h1.fasta ../juicer/aligned/merged_nodups.txt
#after manually curation
run-asm-pipeline-post-review.sh -i 50000 -r Sherbacea_v2.h1.final2.review.assembly ../juicer/reference/Sherbacea_v2.h1.fasta ../juicer/aligned/merged_nodups.txt
