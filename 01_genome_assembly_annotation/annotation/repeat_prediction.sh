### Repeat prediction using RepeatMasker v4.1.5 based on repeat sequences registered in Repbase, species-specific repeat datasets from Populus tricocarpa, de novo repeat libraries built with RepeatModler/2.0.3 83 and LTR_retriever 

# 1- de novo 
module load bioinfo-tools RepeatModeler/2.0.3
RepeatModeler -database ../../index/Sherbacea_final.h2 -LTRStruct -genomeSampleSizeMax 330000000 -engine ncbi -pa 10 &> Sherbacea_final.h2.repeatmodeler.out

# 2- integrate database iteratively 
module load bioinfo-tools bioawk
module load bioinfo-tools RepeatModeler/2.0.3
module load bioinfo-tools RepeatMasker/4.1.5
module load bioinfo-tools SeqKit/2.4.0

cat Sherbacea_final.h2-families.fa | seqkit fx2tab | grep -v "Unknown" |seqkit tab2fx > Sherbacea_final.h2-families.known.fa
cat Sherbacea_final.h2-families.fa | seqkit fx2tab | grep "Unknown" |seqkit tab2fx > Sherbacea_final.h2-families.unknown.fa

/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -t 10 -u Sherbacea_final.h2-families.unknown.fa -k Sherbacea_final.h2-families.known.fa -a Sherbacea_final.h2-families.known.fa -o round1_h2-self
/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -d Arabidopsis -t 10 -u round1_h2-self/round1_h2-self.unknown -k round1_h2-self/round1_h2-self.known -a round1_h2-self/round1_h2-self.known -o round2_h2-self
/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -d Arabidopsis -t 10 -u round2_h2-self/round2_h2-self.unknown -k round2_h2-self/round2_h2-self.known -a round2_h2-self/round2_h2-self.known -o round3_h2-self
/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -d Arabidopsis -t 10 -u round3_h2-self/round3_h2-self.unknown -k round3_h2-self/round3_h2-self.known -a round3_h2-self/round3_h2-self.known -o round4_h2-self
/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -d Arabidopsis -t 10 -u round4_h2-self/round4_h2-self.unknown -k round4_h2-self/round4_h2-self.known -a round4_h2-self/round4_h2-self.known -o round5_h2-self
/proj/snic2021-6-33/bioinformatics_tools/GenomeAnnotation/repclassifier -d Arabidopsis -t 10 -u round5_h2-self/round5_h2-self.unknown -k round5_h2-self/round5_h2-self.known -a round5_h2-self/round5_h2-self.known -o round6_h2-self