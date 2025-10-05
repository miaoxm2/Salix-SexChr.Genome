## 1- prepare alignment matrix

module load bioinfo-tools muscle/5.1
# input file
PAIR_FILE="gene.pair"
PRO_FA="pro.fa"
CDS_FA="cds.fa"
OUTDIR="msa-pal2nal.chr15.h1-h2"

# outdir
mkdir -p "$OUTDIR"

# read each pairs
while read -r line; do
    gene1=$(echo "$line" | awk '{print $1}')
    gene2=$(echo "$line" | awk '{print $2}')
    base_name="${gene1}-${gene2}"
    echo "Processing: $base_name"

    echo $gene1 > tmp.genepair
    echo $gene2 >> tmp.genepair

    # extract seqs
    pullseq -n tmp.genepair -i "$PRO_FA" > tmp.pro.fa
    pullseq -n tmp.genepair -i "$CDS_FA" > tmp.cds.fa

    # align pro seq
    muscle -align tmp.pro.fa -output tmp.pro.aln.fa

    # align cds based on pro
    /proj/snic2021-6-33/bioinformatics_tools/pal2nal.v14/pal2nal.pl tmp.pro.aln.fa tmp.cds.fa -output fasta  -nomismatch -nogap > $OUTDIR/${base_name}.aligned.fasta

done < "$PAIR_FILE"

## 2- build database
genepair="/crex/proj/snic2021-6-33/0xiaomeng/14-geneloss/1-homo-zw/sher.h1chr15-h2chr15.pro.hard.id80.pairs"

pro1="/crex/proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/Sherbacea_final.h1.chr15.pro.fasta"
pro2="/crex/proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/Sherbacea_final.h2.chr15.pro.fasta"

cds1="/crex/proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/Sherbacea_final.h1.chr15.cds.fasta"
cds2="/crex/proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/Sherbacea_final.h2.chr15.cds.fasta"


cd /crex/proj/snic2021-6-33/0xiaomeng/15-deleterious-variation/paml/1-h1h2
cat $pro1 $pro2 > pro.fa
cat $cds1 $cds2 > cds.fa
sed 's/=/\t/g' $genepair > gene.pair
echo "10" > proc

## 3- infer Ks Ka/Ks
# use muscle remove gap/mismatch
ParaAT.pl -h gene.pair -n cds.fa -a pro.fa -p proc -m muscle -f axt -o msa-axt.v2.chr15.h1-h2 -g -t
