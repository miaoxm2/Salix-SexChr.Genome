## 02_sex-linked-region

### **2-1 Genotype calling of females and males**

1. Whole-genome short reads for both pools and individuals were preprocessed with Fastp 
2. We mapped clean paired-end reads to our genome assembly using BWA v0.7.17. 
3. The mapping output was filtered and sorted according to coordinates using SAMtools v1.17. 
4. Duplicate reads from sample preparation, amplification and sequencing were detected and marked using PicardTools. 
5. Variants of individual samples were called and genotyped using GATK v4.1.4.1.
6. We applied hard filtering pipelines to SNPs.
7. For the pooled samples, the variants were called using SAMtools v1.17, and PoPoolation2 were then used to call SNPs.

### **2-2 Identification of sex-linked regions**
We utilized depth-based, genetic differentiation-based and heterozygosity-based approaches.
1. Mapping depth of females and males in 10kb windows was calculated using Mosdepth.
2. Genetic differentiation as weighted Weir and Cockerham’s FST between sexes in 10kb windows with 10kb steps was calculated using VCFtools.
3. The proportion of heterozygous sites in 10kb windows of each re-sequenced individual was calculated manually.
4. Sex-specific regions were identified based on read depth.
5. Female-specific polymorphisms were extracted by comparing female and male genotypes.

### **2-3 Structural variation on chr. 15 in _S. herbacea_**
1. Haplotypes within and between species were aligned using the MUMmer v3.23 pipeline.
2. LTR_harvest and LTR_finder were used to predict full-length LTR retrotransposons in the genome assembly.

### **2-4 Pseudogenization and gene loss**
1. Pseudogenes were annotated using a modified version of the pipelines of Xie et al.
2. Conserved genes, ancestral gene, and conserved/ancestral genes lost on each haplotype

### **2-5 Degeneration and divergence**
1. Variants were annotated using SnpEff v5.2
2. We estimated the synonymous site divergence rates (Ks) and the ratio between nonsynonymous and synonymous divergence rates (Ka/Ks) of homologous genes between the Z- and W-haplotypes based on the YN method using KaKs_caculator 2.0.

