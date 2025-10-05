module load bioinfo-tools LTR_Finder/1.0.7
module load bioinfo-tools LTR_retriever/2.9.0
#predict LTR
fasta="Salix.herbacea/Sherbacea_final.h1.fasta"
LTR_FINDER_parallel -seq $fasta -t 10 -harvest_out -w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85


module load bioinfo-tools RepeatMasker/4.1.5
module load bioinfo-tools blast/2.15.0+
module load bioinfo-tools hmmer/3.4
module load bioinfo-tools trf/4.09.1
module load bioinfo-tools  GenomeTools/1.6.1


gt suffixerator -db Sherbacea_final.h1.fasta -indexname Sherbacea_final.h1.fasta -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Sherbacea_final.h1.fasta -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > Sherbacea_final.h1.fasta.harvest.scn
cat Sherbacea_final.h1.fasta.harvest.scn Sherbacea_final.h1.fasta.finder.combine.scn > Sherbacea_final.h1.fasta.rawLTR.scn

#estimate LTR insertion time
LTR_retriever -genome Sherbacea_final.h1.fasta -inharvest Sherbacea_final.h1.fasta.rawLTR.scn -threads 10

LAI -genome Sherbacea_final.h1.fasta -intact Sherbacea_final.h1.fasta.pass.list -all Sherbacea_final.h1.fasta.out
