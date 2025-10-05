
### prepare files. extract seq id from single copy orthologues
bash ~/scripts/extract_genes_from_gff.sh Orthogroups_SingleCopyOrthologues.txt Orthogroups.txt singlecopy.orthogroups
grep -v -e 'chr15' -e 'chr07' -e 'chr19' singlecopy.orthogroups > singlecopy.orthogroups.rmchr07_15_19

### extract fasta sequence.per orthogroups
# list the fasta name according to the order in orthogroups.
### create MSA using coding sequences

pqio="/0-Salicaceae_genome/Populus.qiongdaoensis/Pqiongdaoensis.cds.fasta"
ptri="/0-Salicaceae_genome/Populus.trichocarpa/ptri.cds.fa"
sarb=/0-Salicaceae_genome/Salix.arbutifolia/sarb_a.cds.fa
sbra="/0-Salicaceae_genome/Salix.brachista/sbra.cds.fa"
sdun="/0-Salicaceae_genome/Salix.dunnii/sdun.cds.fa"
sher="/0-Salicaceae_genome/Salix.herbacea/sher_h1.cds.fa"
skor="/0-Salicaceae_genome/Salix.koriyanagi/skor.cds.fa"
spur=/0-Salicaceae_genome/Salix.purpurea/1-femaleclone-94006/spur.cds.fa
ssuc=/0-Salicaceae_genome/Salix.suchowensis/ssuc_xy12f.cds.fa
sude="/0-Salicaceae_genome/Salix.udensis/sude.cds.fa"
svim="/0-Salicaceae_genome/Salix.viminalis/svim_map2.cds.fa"

cat $pqio $ptri $sarb $sbra $sdun $sher $skor $spur $ssuc $sude $svim > 11sp.cds.fa
cat ../../../*fa > 11sp.pro.fa

awk '{$1=""; sub(/^ /, ""); print}' ../Orthogroups/singlecopy.orthogroups.rmchr07_15_19 |sed 's/ /\t/g' > 11sp.sco.genelist

module load bioinfo-tools  MAFFT/7.407
echo "10" > proc
rm -rf msa-cds
ParaAT.pl -h 11sp.sco.genelist -n 11sp.cds.fa -a 11sp.pro.fa -o msa-cds -p proc -msa mafft -format paml -nogap -nomismatch

### transfer alignment output into paml-readable phylip

perl merge_phylip.split_codon.pl 'msa-cds/*paml' 11sp.aligned_cds.phy


### manually check if any error infos. and then cat them
ls *aligned.rename.fa.treefile |wc -l
ls *aligned.rename.fa | wc -l

cat msa/*aligned.rename.fa.treefile > sco_1115tree.txt

astral="/crex/proj/snic2021-6-33/bioinformatics_tools/Astral/astral.5.7.8.jar"
java -jar $astral -i sco_1115tree.txt -o sco_1115tree.species_tree.txt 2>out.log

### run mcmctree
module load bioinfo-tools paml/4.9j
mcmctree mcmctree.ctl

