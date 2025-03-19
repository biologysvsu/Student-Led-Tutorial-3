# Actual tutorial
## Download BEFORE Class
### Download R and RStudio
  - [R Project](https://www.r-project.org/)

#### Steps to download: Under Download header > CRAN > 0-Cloud > download for your computer type (should be a .pkg file) > install > open R and R Studio

# Student-Led-Tutorial-3
# Task: Tutorial for Gene Expression Analysis Using STAR, Gene Expression Quantification, and Differential Expression Analysis
# Date: March 20th

## **Objective**
Students will:
1. Download RNA-seq data for **mock-infected** and **COVID-19-infected** human cell lines from the provided SRA project.
2. Perform RNA-seq alignment using **STAR**.
3. Quantify gene expression using **FeatureCounts**.
4. Perform differential expression analysis between mock-infected and COVID-19-infected cells using **DESeq2** (or equivalent).
5. Visualize results such as heatmaps and volcano plots to highlight key differentially expressed genes.

---
# Load Modules

  ``` bash
  module load anaconda3/2024.10-1
  conda activate bioinfo-env
  module load STAR 

