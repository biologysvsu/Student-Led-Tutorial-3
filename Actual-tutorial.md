# Actual tutorial
## Download BEFORE Class
### Download R and RStudio
  - [R Project](https://www.r-project.org/)

#### Steps to download: Under Download header > CRAN > 0-Cloud > download for your computer type (should be a .pkg file) > install > open R and R Studio

# Student-Led-Tutorial-3

# Task: Tutorial for Gene Expression Analysis Using STAR, Gene Expression Quantification, and Differential Expression Analysis
  ## In simpler terms: We are going to use RNA-seq analysis for reading single-end reads and then visualize the results
  
# Date: March 20th

## **Objective**
Students will:
1. Download RNA-seq data for **mock-infected** and **COVID-19-infected** human cell lines from the provided SRA project.
2. Perform RNA-seq alignment using **STAR**.
3. Quantify gene expression using **FeatureCounts**.
4. Perform differential expression analysis between mock-infected and COVID-19-infected cells using **DESeq2** (or equivalent).
5. Visualize results such as heatmaps and volcano plots to highlight key differentially expressed genes.
---
**Part 1: Setup & Environment Preparation**
## Step 1: Navigate to the shared directory and create a new  working directory for this tutorial:
  ```bash
  cd /ocean/projects/agr250001p/shared/ 
  mkdir Tutorial-3 
  cd Tutorial-3
  ```
## Step 2: Load anaconda and create bioinformatics environment
```bash
module load anaconda3/2024.10-1
conda create --name bioinfo-env python=3.9 -y 
conda activate bioinfo-env 
conda install -c bioconda star subread fastqc samtools -y
```

**Part 2: Data Prepatration**
## Step 1: Download FASTQ Files
  ### We already extracted these for you. If needed, verify:
  ```bash
  ls -lh *.fastq
  ```
## Step 2: Download Necessary Files
## Download the Reference Genome
  ```bash
  wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
  ```
### Gunzip the 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' file
  ```bash
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
  ```
## Download the Annotation File 
  ```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
  ```
### Gunzip the 'Homo_sapiens.GRCh38.113.gtf.gz' file
  ```bash
gunzip Homo_sapiens.GRCh38.113.gtf.gz
  ```
## Move the 'Homo_sapiens.GRCh38.113.gtf' file into the 'annotation.gtf' file
 ```bash
mv Homo_sapiens.GRCh38.113.gtf annotation.gtf
  ```

STAR --runMode genomeGenerate --genomeDir /path/to/genome --genomeFastaFiles /path/to/genome.fasta --sjdbGTFfile /path/to/annotation.gtf --runThreadN 8
```
### You now should have a STAR index
- To check:
  ```bash
  ls -lh STAR_index
  ```
## Create STAR Alignment Script
  ```bash
  nano star_alignment.sbatch
  ```
### Paste the following content into "star_alignment.sbatch':
  ```bash
#!/bin/bash #SBATCH --job-name=STAR_alignment #SBATCH --output=star_alignment.out #SBATCH --error=star_alignment.err #SBATCH --mem=50G #SBATCH --time=4:00:00 #SBATCH --cpus-per-task=16 #SBATCH --mail-type=END,FAIL #SBATCH --mail-user=jrander3@svsu.edu # Change to your email 

module load STAR 

for file in SRR11412215.fastq.gz SRR11412216.fastq.gz SRR11412217.fastq.gz SRR11412218.fastq.gz  

SRR11412227.fastq.gz SRR11412228.fastq.gz SRR11412229.fastq.gz SRR11412230.fastq.gz SRR11412231.fastq.gz do STAR --genomeDir star_index
  ```
### You have now completed a STAR alignment for the Mock infected cells and
  
