### scaffolding using HiC reads

module load bioinfo-tools lastz/1.04.00
#rm -f Sherbacea_v2.h1.*.cprops Sherbacea_v2.h1.mnd.*.txt
run-asm-pipeline.sh -i 50000 ../juicer/reference/Sherbacea_v2.h1.fasta ../juicer/aligned/merged_nodups.txt

#use JuiceBox to manually curate
#after manually curation
module load bioinfo-tools lastz/1.04.00
run-asm-pipeline-post-review.sh -i 50000 -r Sherbacea_v2.h1.0.review_finalize3.assembly ../0-juicer/reference/Sherbacea_v2.h1.fasta ../0-juicer/aligned/merged_nodups.txt