# WWWG2B Project
All codes associated with the WWWG2B project are organized in this repository. Below is a guide to each directory and the codes it contains:

## Directory Structure
```
WWWG2B
|-- 00.SNP_calling_and_QC
|-- 01.HAPMAP_pipeline
|-- 02.PCA_and_t-SNE
|-- 03.NAM_imputation_pipeline
|-- 04.GWAS_pipeline
|-- 05.HAPPE_for_wheat
|-- 06.tagSNP_pipeline
|-- 07.CNV_pipeline
```

---

### 00.SNP_calling_and_QC
This directory contains scripts for performing the entire SNP calling workflow, starting from BWA alignment to SNP calling. Additionally, it contains filtering scripts for SNPs and INDELs.

#### Key Features
- BWA alignment
- SNP calling
- SNP & INDEL filtering

---

### 01.HAPMAP_pipeline
This directory contains code for creating HAPMAP files, specifically scripts to connect blocks during the construction phase of HAPMAP.

#### Key Features
- Block connecting
- HAPMAP construction

---

### 02.PCA_and_t-SNE
In this directory, you will find code for Principal Component Analysis (PCA) and t-Distributed Stochastic Neighbor Embedding (t-SNE) of SNP data.

#### Key Features
- PCA on SNP data
- t-SNE on SNP data

---

### 03.NAM_imputation_pipeline
This directory contains scripts for imputing data in NAM populations. It includes two versions: a single-threaded version and a multi-threaded version.

#### Key Features
- Single-threaded imputation
- Multi-threaded imputation

---

### 04.GWAS_pipeline
This directory contains scripts for running Genome-Wide Association Studies (GWAS) analyses. The pipeline supports multiple software options and includes visualization scripts for Manhattan and QQ plots.

#### Key Features
- GWAS analysis
- Manhattan plot generation
- QQ plot generation

---

### 05.HAPPE_for_wheat
This directory is specifically designed for wheat genomes and contains a specialized version of HAPPE. For a general-purpose version of HAPPE, please refer to [this GitHub repository](https://github.com/fengcong3/HAPPE). Additionally, it contains R scripts for visualizing haplotype and depth data.

#### Key Features
- Specialized HAPPE for wheat
- Haplotype visualization
- Depth data visualization

---
### 06.tagSNP_pipeline
This directory contains scripts for selecting tagSNPs from a SNP matrix using HAPMAP data. These selected, minimal SNPs are capable of distinguishing between 1,000 samples and can be utilized for marker design.

#### Key Features
- TagSNP selection from SNP matrix
- Utilizes HAPMAP data
- Suitable for marker design in large datasets

---

---
### 07.CNV_pipeline
This directory contains scripts for calling Copy Number Variants (CNVs) from bam file.

#### Key Features
- CNV calling

---

Feel free to dive into each directory for more details on the individual codes and pipelines.
