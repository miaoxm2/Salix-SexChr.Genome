### gene prediction following automated pipelines using BRAKER with GeneMark-ETP, AUGUSTUS and TSEBRA for three separate rounds: homology-based only, RNA-seq only, and combining homology-based and RNA-seq

#1-rna-align/mapping-hisat.sh
r1="/proj/snic2021-6-33/0xiaomeng/9-annotation/2-RNA/rna_1_clean.fq.gz"
r2="/proj/snic2021-6-33/0xiaomeng/9-annotation/2-RNA/rna_2_clean.fq.gz"
ref="/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/index/Sherbacea_final.h1"
out="rna-h1_masked"
module load bioinfo-tools HISAT2/2.2.1
module load bioinfo-tools samtools/1.17
module load bioinfo-tools picard/2.27.5

##build library in index/
hisat2-build Sherbacea_final.h1.full_mask.soft.fasta Sherbacea_final.h1
##mapping
hisat2 -p 10 -5 3 -3 3 --mp 3,0 -x $ref -1 $r1 -2 $r2 -S $out.sam 
samtools view -f2 -q20 -b -@10 $out.sam | samtools sort -@10 - > $out.bam
java -jar $PICARD_ROOT/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=$ref INPUT=$out.bam OUTPUT=$out.bam_alignment.stats


#2-braker-rna/braker.rna.sh
module load bioinfo-tools braker/2.1.6
source $AUGUSTUS_CONFIG_COPY
#cp -vf /sw/bioinfo/GeneMark/keyfile/gm_key $HOME/.gm_key
braker.pl --cores 12 --workingdir=try1/ \
 --genome=/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/index/Sherbacea_final.h1.full_mask.soft.fasta  \
 --softmasking --bam=/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/1-rna-align/rna-h1_masked.bam \
 --AUGUSTUS_CONFIG_PATH=/home/maoxm/medium_storage/snic2021-6-33/0xiaomeng/9-annotation/00-final/augustus_config \
 -gff3 --makehub --species=sherbacea --email=xiaomeng.mao@ebc.uu.se

#3-braker-pro/braker.pro.sh
module load bioinfo-tools braker/2.1.6
source $AUGUSTUS_CONFIG_COPY
cp -vf /sw/bioinfo/GeneMark/keyfile/gm_key $HOME/.gm_key
braker.pl --cores 12 --workingdir=try1/ \
  --genome=/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/index/Sherbacea_final.h1.full_mask.soft.fasta \
  --softmasking --prot_seq=/proj/snic2021-6-33/0xiaomeng/9-annotation/3-homology/all_pep.fa \
  --AUGUSTUS_CONFIG_PATH=/home/maoxm/medium_storage/snic2021-6-33/0xiaomeng/9-annotation/00-final/3-braker-pro/augustus_config \
  --species=sherbacea -gff3

#4-braker-rna-pro/braker.prorna.sh
module load bioinfo-tools braker/2.1.6
source $AUGUSTUS_CONFIG_COPY
braker.pl --cores 12 --workingdir=try1/ \
 --genome=/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/index/Sherbacea_final.h1.full_mask.soft.fasta  \
 --prot_seq=/proj/snic2021-6-33/0xiaomeng/9-annotation/3-homology/all_pep.fa \
 --softmasking --bam=/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/1-rna-align/rna-h1_masked.bam \
 --AUGUSTUS_CONFIG_PATH=/home/maoxm/medium_storage/snic2021-6-33/0xiaomeng/9-annotation/00-final/4-braker-rna-pro/augustus_config \
  --etpmode -gff3 --species=sherbacea 

#5-combine/combine.sh
rna="/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/2-braker-rna/try1"
pro="/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/3-braker-pro/try1"
rnapro="/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/4-braker-rna-pro/try1"
out="h1_final.braker_combined.try1"

module load bioinfo-tools TSEBRA/1.0.3-20211126-336c380
#tsebra.py -g $rna/augustus.hints.gtf,$pro/augustus.hints.gtf,$rnapro/braker.gtf  -c $TSEBRA_ROOT/config/default.cfg \
    -e $rna/hintsfile.gff,$pro/hintsfile.gff,$rnapro/hintsfile.gff \
    -o $out.gtf

module load bioinfo-tools TSEBRA/1.0.3-20211126-336c380
rename_gtf.py --gtf $out.gtf --translation_tab translation.tab --out $out.rename.gtf
fix_gtf_ids.py --gtf  $out.rename.gtf --out $out.renamed.fixed.gtf

ref="/proj/snic2021-6-33/0xiaomeng/9-annotation/00-final/index/Sherbacea_final.h1.full_mask.soft.fasta"
module load bioinfo-tools gffread/0.12.6
gffread $out.renamed.fixed.gtf -g $ref -y $out.pro.fa -x $out.cds.fa

module load bioinfo-tools AGAT/1.0.0
agat_convert_sp_gxf2gxf.pl --gtf $out.renamed.fixed.gtf -o $out.renamed.fixed.gff


