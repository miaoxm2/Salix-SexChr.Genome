### evolutionary history of sex-linked genes among Salicaceae
module load bioinfo-tools last
### prepare for $sp data
cds="path/to/codingseq"
sp="species"

python -m jcvi.formats.fasta format $cds $sp.cds
sed -i 's/ Protein.*//g' $sp.cds
python -m jcvi.formats.gff bed --type=mRNA --key=Accession /crex/proj/snic2021-6-33/private/genome/Populus.qiongdaoensis/GWHBJCS00000000.gff.gz -o $sp.bed

python -m jcvi.compara.catalog ortholog $sp sher --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple $sp.sher.anchors $sp.sher.anchors.new

python -m jcvi.compara.synteny mcscan $sp.bed $sp.sher.anchors.new --iter=1 -o $sp.sher.i1.blocks

python -m jcvi.formats.bed merge sher.bed $sp.bed -o sher.$sp.bed

awk 'FNR==NR { gene[$1]++; next } ($2 ~ "chr15"){OFS="\t";if($2 in gene){print "g*"$1,$2} else {print $1,$2}}' list $sp.sher.i1.blocks > $sp.sher.15only.blocks

python -m jcvi.graphics.synteny $sp.sher.15only.blocks $sp.sher.bed $sp.sher.chr15.layout --glyphcolor=orthogroup --notex



### create three species panels

#spur-sher-sbra
python -m jcvi.graphics.synteny spur.sher.15only.blocks spur.sher.bed spur.sher.chr15.layout --glyphcolor=orthogroup --notex

echo "# x,   y, rotation,   ha,     va,   color, ratio, label
0.5, 0.6,        0, left, center, #fc8d62,     1, sher chr15w
0.5, 0.4,        0, left, center,       m,     1, sbra chr15w
# edges
e, 0, 1" > sher.sbra.chr15.layout 
python -m jcvi.graphics.synteny sbra.sher.15only.blocks sbra.sher.bed sher.sbra.chr15.layout --glyphcolor=orthogroup --notex

#sarb-sher-sdun
python -m jcvi.graphics.synteny sarb.sher.15only.blocks sarb.sher.bed sarb.sher.chr15.layout --glyphcolor=orthogroup --notex
python -m jcvi.graphics.synteny sdun.sher.15only.blocks sdun.sher.bed sher.sdun.chr15.layout --glyphcolor=orthogroup --notex

#ptri-sher-pqio
echo "# x,   y, rotation,   ha,     va,   color, ratio, label
 0.5, 0.6,        0, left, center, #fc8d62,     1, sher chr15w
 0.5, 0.4,        0, left, center,       m,     1, pqio chr15
 # edges
 e, 0, 1" > sher.pqio.chr15.layout 

python -m jcvi.graphics.synteny ptri.sher.15only.blocks ptri.sher.bed ptri.sher.chr15.layout --glyphcolor=orthogroup --notex
python -m jcvi.graphics.synteny pqio.sher.15only.blocks pqio.sher.bed sher.pqio.chr15.layout --glyphcolor=orthogroup --notex
