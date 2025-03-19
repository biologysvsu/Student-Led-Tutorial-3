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

**Part 3: STAR Genome Indexing**

## Step 1: Create STAR Indexing Script
  ### Run:
  ```bash
nano star_indexing.sbatch 
  ```

  ### Paste:
  ```bash
#!/bin/bash 
#SBATCH --job-name=star_indexing 
#SBATCH --output=star_indexing.out 
#SBATCH --error=star_indexing.err 
#SBATCH --time=4:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=100G 
#SBATCH --mail-user=jrander3@svsu.edu 
#SBATCH --mail-type=ALL 
  ```

## Step 2: Load STAR
  ```bash
module load STAR 
  ```
  ### Run:
  ```bash
  STAR --runMode genomeGenerate \ 
     --genomeDir star_index \ 
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \ 
     --sjdbGTFfile annotation.gtf \ 
     --sjdbOverhang 100 \ 
     --runThreadN 16

### You now should have a STAR index
- To check:
  ```bash
  ls -lh STAR_index
  ```

## Step 3: Submit STAR Indexing Job
```bash
sbatch star_indexing.sbatch
 ```

## Step 4: Verify Indexing Output
  ```bash
## Create STAR Alignment Script
  ```bash
  nano star_alignment.sbatch
  ```

## Step 5: Create Symlink to Shared STAR Index
- If STAR indexing has already been generated in a shared directory, create a symlink instead:  
  ```bash
  ln -s /ocean/projects/agr250001p/shared/tutorial-data/tutorial_7_data/star_index . 
  ```
  
**Part 4: Read Alignment with STAR**

## Step 1: Align reads to the Genome
  ### For each sample (mock or COVID-infected), align reads to the genome:
  ```bash
  STAR --genomeDir star_index \ 
     --readFilesIn SRR11412215.fastq \ 
     --outFileNamePrefix mock_rep1_ \ 
     --outSAMtype BAM SortedByCoordinate  
  ```
  ### Replace 'SRR11412215' with the appropriate sample name for each replicate
  
## Step 2: Expected Output Files
  ### Each Sample should generate:
  - Sorted BAM file: mock_rep1_Aligned.sortedByCoord.out.bam 
  - Alignment log file: mock_rep1_Log.final.out
   
## Step 3: Setup and Load Required Libraries
### Ensure that all required packages are installed and loaded:
```bash
# Install Bioconductor and DESeq2 if not installed 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
BiocManager::install("DESeq2") 
 
# Load required libraries 
library(DESeq2) 
library(ggplot2) 
library(pheatmap) 
  ```
## Step 4: Load and Prepared Count Data
  - Verify File location
  - Set the working directory or use file.choose() to select the count file:
```bash
    # Set working directory if needed 
setwd("C:/Users/19897/Desktop/") 
 
# Verify file exists 
file.exists("counts.txt")  # Should return TRUE 
```
## Step 5: Read Count Data
## Load count data from the file:
  ```bash
counts <- read.table("C:/Users/19897/Desktop/counts.txt", 
                     header=TRUE, 
                     row.names=1, 
                     sep="\t", 
                     check.names=FALSE)
  ```
## Remove metadata Columns (Chr, Start, End, Strand, Length) 
  ```bash
counts <- counts[, -(1:5)]
  ```
## Verify data structure
- 'dim(counts)' checks dimensions (genes, samples)
- 'colnames(counts)' ensures only sample names remain
```bash
dim(counts)   
```
```bash
colnames(counts)   
```
**Part 5: Create Sample Metadata**
## Define experimental conditions for each sample: 
  ```bash
coldata <- data.frame( 
  condition = c(rep("mock", 4), rep("covid", 5)),  # Adjust sample numbers if needed 
  replicate = c(1:4, 1:5)  )  
  ```
## Assign sample names
```bash
rownames(coldata) <- colnames(counts)
```
## Check that row names match column names (Response should say: TRUE)
```bash
all(rownames(coldata) == colnames(counts))
```
**Part 6: Differential Expression Analysis with DESeq2**
## Create DESeq2 dataset
```bash
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
```
## Run DESeq2 differential expression analysis
```bash
dds <- DESeq(dds)
```
## Get results
```bash
results <- results(dds
```
## Save results to CSV
```bash
write.csv(as.data.frame(results), "differential_expression_results.csv")
```
**Part 7: Visualize Results**
## Step 1: Volcano Plot
  - Identify Significant genes
```bash
results$significant <- !is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1
```
## Step 2: Generate Volcano Plot
```bash
ggplot(results, aes(x=log2FoldChange, y=-log10(padj), color=significant)) + 
  geom_point(alpha=0.7) + 
  theme_minimal() + 
  scale_color_manual(values=c("gray", "red")) + 
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10 Adjusted p-value")
```
## Step 3: Select the top 50 genes based on adjusted p-value
```bash
library(pheatmap) 
top_genes <- rownames(results)[order(results$padj, na.last=NA)][1:50]
```
## Step 4: Normalize couns for heatmap visualization
```bash
normalized_counts <- counts(dds, normalized=TRUE)
```
## Step 5: Generate heatmap
```bash
pheatmap(normalized_counts[top_genes, ], cluster_rows=TRUE, cluster_cols=TRUE, scale="row")
```

#**The End**
