### scaffolding using HiC reads

wkpath='/proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/0-juicer'
#software= "/proj/snic2021-6-33/bioinformatics_tools/juicer-1.6/juicer" & "juicer_modifiled.sh"

REF='/proj/snic2021-6-33/0xiaomeng/5-call_variants/00-version2/index/Sherbacea_v2.h1.fasta'
assembly='she_h1'
enzyme='HindIII'
FASTQ1='/proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/clean_reads/BC425-01H0001_1.filter3.fastq'
FASTQ2='/proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/clean_reads/BC425-01H0001_2.filter3.fastq'
SRA='BC425-01H0001'

REP_LABEL=$assembly\_rep2

module load bioinfo-tools bwa/0.7.17

echo 'step0-1.prepare-index'
mkdir references/ && cd references/
ln -s $REF $assembly.fasta
bwa index $assembly.fasta
cd ..

echo 'step0-2.prepare-fastq'
mkdir fastq/ 
cd fastq/
ln -s $FASTQ1 $SRA\_R1.fastq
ln -s $FASTQ2 $SRA\_R2.fastq
cd ..
###this step will creat too many splited files
mkdir splits/
cd splits/
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/$SRA\_R2.fastq &
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/$SRA\_R1.fastq &
wait
cd ..

echo 'step0-3.prepare-emzyme'
mkdir restriction_sites/
cd restriction_sites/
python /proj/snic2021-6-33/bioinformatics_tools/juicer-1.6/misc/generate_site_positions.py $enzyme $assembly ../references/$assembly.fasta 
awk 'BEGIN{OFS="\t"}{print $1, $NF}' $assembly\_$enzyme.txt > $assembly.chrom.sizes
cd ..

echo 'step0-4.prepare-scripts'
ln -s /proj/snic2021-6-33/bioinformatics_tools/juicer-1.6/CPU scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
cd ../..

echo 'step1.juicer'
rm -rf aligned/
./juicer_modifiled-D.sh -g $REP_LABEL -z references/$assembly.fasta \
-p restriction_sites/$assembly.chrom.sizes \
-y restriction_sites/$assembly\_$enzyme.txt\
-s HindIII \
-d /proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/juicer2 \
-t 10 > log.txt
