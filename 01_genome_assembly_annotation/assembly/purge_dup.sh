### purge duplicates
#Step 1-1. map reads to assembly
module load bioinfo-tools minimap2/2.24-r1122

minimap2 -xasm20 -t 12 /proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/hifiasm/run2/hifi_hic_v2_s3.hap1.fasta /proj/snic2021-6-33/0xiaomeng/0-Salicaceae_genome/Salix.herbacea/fastq_F005/pt_056_001_all.filter.fastq | gzip -c - > she.hifiasm.purge.paf.gz
/proj/snic2021-6-33/bioinformatics_tools/purge_dups/bin/pbcstat she.hifiasm.purge.paf.gz 
/proj/snic2021-6-33/bioinformatics_tools/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log

#Step 1-2. map reads to assembly

/proj/snic2021-6-33/bioinformatics_tools/purge_dups/bin/split_fa /proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/hifiasm/run2/hifi_hic_v2_s3.hap1.fasta > she.hifi_hic_h1.split
minimap2 -xasm5 -t 12 -DP she.hifi_hic_h1.split she.hifi_hic_h1.split | gzip -c - > she.hifi_hic_h1.split.self.paf.gz

#Step 2. Purge haplotigs and overlaps with the following command.
/proj/snic2021-6-33/bioinformatics_tools/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov she.hifi_hic_h1.split.self.paf.gz > dups.bed 2> purge_dups.log

#Step 3. Get purged primary and haplotig sequences from draft assembly.
/proj/snic2021-6-33/bioinformatics_tools/purge_dups/bin/get_seqs -e dups.bed /proj/snic2021-6-33/0xiaomeng/10-assembly/1hic/hifiasm/run2/hifi_hic_v2_s3.hap1.fasta 
