## construct alignment among genes
module load bioinfo-tools clustalw/2.1
echo "10" > proc
ParaAT.pl -h ../Orthogroups/singlecopy.genepairs -n ptri-sher-spur.cds.fasta -a ptri-sher-spur.chr15.pro.fasta -p proc -f fasta -o msa.ptri-sher-spur -m clustalw2

## build phylogenetic tree
module load bioinfo-tools iqtree
for fasta in msa.ptri-sher-spur/*fasta; do
iqtree2 -s $fasta -m MFP -st CODON -T AUTO
done

## cluster trees in R
