## detect ortholog among ptri.chr15.pro.fasta; sarb.chr15y.pro.fasta; sher.chr15w.pro.fasta; spur.chr15w.pro.fasta; ssuc.chr15w.pro.fasta; pqio.chr15.pro.fasta; sarb.chr15x.pro.fasta; sdun.chr15.pro.fasta; sher.chr15z.pro.fasta; spur.chr15z.pro.fasta

module load bioinfo-tools OrthoFinder/2.5.2
orthofinder.py -t 15 -f 4-versionfinal-chrz/

## find h1/h2 orthologous with any of other species
awk '($7>0 && $7!=$12 && $7+$8<$12){print $1}' Orthogroups.GeneCount.tsv > h1-shared_w_anyother.orth
awk '($8>0 && $8!=$12 && $7+$8<$12){print $1}' Orthogroups.GeneCount.tsv > h2-shared_w_anyother.orth

bash ~/scripts/extract_genes_from_gff.sh h1-shared_w_anyother.orth Orthogroups.txt h1-shared_w_anyother.orthgroups
bash ~/scripts/extract_genes_from_gff.sh h2-shared_w_anyother.orth Orthogroups.txt h2-shared_w_anyother.orthgroups

#check
wc -l h1-shared_w_anyother.orth*
wc -l h2-shared_w_anyother.orth*
#generate geneid
awk '{ for (i=1; i<=NF; i++) if ($i ~ /Saher.h1_chr15/) print $i }' h1-shared_w_anyother.orthgroups | sort | uniq > h1-shared_w_anyother.orthgroups.gene
awk '{ for (i=1; i<=NF; i++) if ($i ~ /Saher.h2_chr15/) print $i }' h2-shared_w_anyother.orthgroups | sort | uniq > h2-shared_w_anyother.orthgroups.gene


## find h1/h2 orthologous with all of other three clades of species : Populus/XYspecies/ZWspecies
awk '(($2>0 || $3>0) && ($4>0 || $5>0 || $6>0) && ($9>0 || $10>0 || $11>0) && $7>0 ){print $1} ' Orthogroups.GeneCount.tsv > h1-shared_w_all.orth
bash ~/scripts/extract_genes_from_gff.sh h1-shared_w_all.orth Orthogroups.txt h1-shared_w_all.orthgroups
awk '{ for (i=1; i<=NF; i++) if ($i ~ /Saher.h1_chr15/) print $i }' h1-shared_w_all.orthgroups | sort | uniq > h1-shared_w_all.orthgroups.gene

awk '(($2>0 || $3>0) && ($4>0 || $5>0 || $6>0) && ($9>0 || $10>0 || $11>0) && $8>0 ){print $1} ' Orthogroups.GeneCount.tsv > h2-shared_w_all.orth
bash ~/scripts/extract_genes_from_gff.sh h2-shared_w_all.orth Orthogroups.txt h2-shared_w_all.orthgroups
awk '{ for (i=1; i<=NF; i++) if ($i ~ /Saher.h2_chr15/) print $i }' h2-shared_w_all.orthgroups | sort | uniq > h2-shared_w_all.orthgroups.gene
