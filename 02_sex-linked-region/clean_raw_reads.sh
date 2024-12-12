## QC. clean raw reads
module load bioinfo-tools fastp


for f1 in /proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/fastq_ind/*R1.clean.fastq
do
    f2=${f1%%R1.clean.fastq}"R2.clean.fastq"
    out=/proj/naiss2023-6-57/1-trim_data/
    ff1=${f1##*/}
    ff2=${f2##*/}
    out1=${out}${ff1%%clean.fastq}"clean2.fastq"
    out2=${out}${ff2%%clean.fastq}"clean2.fastq"
    fh1=${out}${ff1%%R1.clean.fastq}"clean.html"

    fastp -w 10 \
    -i $f1 -I $f2 \
    -o $out1 -O $out2 \
    -h $fh1 \
    -q 20 -u 20 -l 100 -y 30 -g 10 -x 10 -p 20

done
