### clean raw data for HiC reads

module load bioinfo-tools fastp/0.23.1
fastp -w 15 -q 20  -i /proj/snic2021-6-33/rawdata/20221229_hic/new/rawdata/Unknown_BC425-01H0001_1.fq.gz  -o BC425-01H0001_1.filter3.fastq \
-I /proj/snic2021-6-33/rawdata/20221229_hic/new/rawdata/Unknown_BC425-01H0001_2.fq.gz -O BC425-01H0001_2.filter3.fastq \
-h BC425-01H00013.html -j BC425-01H00013.json \
--trim_poly_g --poly_g_min_len 10 \
--trim_poly_x --poly_x_min_len 10 \
--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
--qualified_quality_phred 20 \
--length_required 100 \
--detect_adapter_for_pe --correction --dont_overwrite

module load bioinfo-tools FastQC/0.11.9
fastqc -t 15 -o ../fastqc/ BC425-01H0001_1.filter3.fastq BC425-01H0001_2.filter3.fastq 

