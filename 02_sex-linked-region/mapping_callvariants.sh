#! /bin/bash 
ref="/proj/naiss2023-6-57/1-genome_index/Sherbacea_final.h1.softmask.fasta"
for read1 in /crex/proj/naiss2023-6-57/1-trim_seq/bgi2024/X*R1.clean.fastq
do
        read2=${read1/R1/R2}
        file=${read1##/crex/proj/naiss2023-6-57/1-trim_seq/bgi2024/}
        sp=${file%%.R*}
        out="/proj/naiss2023-6-57/2-map_bam/h1mask_She1-"${sp}
        outdup="/proj/naiss2023-6-57/3-markdup/h1mask_She1-"${sp}
        outvcf="/proj/naiss2023-6-57/4-vcf/1-raw_vcf/h1mask_She1-"${sp}
        echo "#! /bin/bash -l 
#SBATCH -A naiss2024-5-173
#SBATCH -p core 
#SBATCH -n 10
#SBATCH -t 2-00:00:00 
#SBATCH -J map.h1She1-${sp}
#SBATCH --mail-user xiaomeng.mao@ebc.uu.se
#SBATCH --mail-type=FAIL,END

module load bioinfo-tools bwa/0.7.17
module load bioinfo-tools samtools/1.17
module load bioinfo-tools picard/2.27.5
###1.mapping, paired
bwa mem -t 10 $ref $read1 $read2 2>>$out.err | samtools view -f2 -bt -@10 -T $ref - -o $out.bam 
###2.sort
samtools sort -@ 10 $out.bam -o $out.sort.bam 2>>$out.sort.log
###3.stats
samtools stats $out.sort.bam > $out.sort.stats
###4.markdup
java -jar \$PICARD_ROOT/picard.jar MarkDuplicates I=$out.sort.bam O=$outdup.markdup.sort.bam M=$outdup.markdup.metrics.txt CREATE_INDEX=True 
java -jar \$PICARD_ROOT/picard.jar AddOrReplaceReadGroups I=$outdup.markdup.sort.bam O=$outdup.addrg.markdup.sort.bam RGID=$sp RGLB=$sp RGPL=illumina RGPU=unit1 RGSM=$sp 
###make index
samtools index -@10 $outdup.addrg.markdup.sort.bam

###5.callsnp
echo '#! /bin/bash -l 
#SBATCH -A naiss2024-5-173
#SBATCH -p core 
#SBATCH -n 4
#SBATCH -t 8-00:00:00 
#SBATCH -J call.h1She1-$sp
#SBATCH --mail-user xiaomeng.mao@ebc.uu.se
#SBATCH --mail-type=END,FAIL
module load bioinfo-tools GATK/4.3.0.0
gatk HaplotypeCaller \
-R $ref \
-I $outdup.addrg.markdup.sort.bam \
-O $outvcf.g.vcf.gz \
-ERC GVCF 2>>$outvcf.gatk.error 
' > 4-callsnp/4-callsnp.h1She1-$sp.sh

" > 2-bwa/script.h1She1-${sp}.sh
sbatch 2-bwa/script.h1She1-${sp}.sh
done

