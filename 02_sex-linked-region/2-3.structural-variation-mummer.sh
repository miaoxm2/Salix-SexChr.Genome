### detect structural variation between haplotypes
ref='sher.h2.chr15.fasta'
query='sher.h1.chr15.fasta'
prefix='h1-h2.chr15.mum'

module load bioinfo-tools MUMmer/3.23
nucmer -g 500 -c 100 -p $prefix $query $ref
delta-filter -i 89 -l 1000 -1 $prefix.delta > $prefix.uniq.delta
show-coords -r $prefix.uniq.delta > $prefix.uniq.coords
mummerplot --fat --color --postscript -p $prefix.uniq -f $prefix.uniq.delta 
ps2pdf $prefix.uniq.ps pdf/$prefix.uniq.pdf
