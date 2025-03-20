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

# PLEASE DO NOT START THE TUTORIAL, DUE TO ISSUES WITH FILE SHARING YOU WILL WATCH THE FIRST 4 PARTS

# Part 1: Setup & Environment Preparation


### Step 1: Navigate to the shared directory and create a new  working directory for this tutorial:
  ```bash
  cd /ocean/projects/agr250001p/shared/ 
  mkdir Tutorial-3 
  cd Tutorial-3

  ```
# Create symlink to the SRA data
```
ln -s /ocean/projects/agr250001p/shared/tutorial-data/tutorial_3_data/sra_data .
```

# Create symlink to the STAR index
```
ln -s /ocean/projects/agr250001p/shared/tutorial-data/tutorial_3_data/star_index .
```

### Step 2: Load anaconda and create bioinformatics environment
```bash
module load anaconda3/2024.10-1
conda create --name bioinfo-env python=3.9 -y 
conda activate bioinfo-env 
conda install -c bioconda star subread fastqc samtools -y
```
Deactivate conda
```
conda deactivate
```

# Part 2: Data Preparation

### Step 1: Download FASTQ Files

  #### We already extracted these for you. If needed, verify:
  ```bash
  ls -lh sra_data
  ```

## Step 2: Download Necessary Files

### Download the Reference Genome
  ```bash
  wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
  ```

#### Gunzip the `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` file
  ```bash
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  ```
#### Gunzip all fastq files
``` bash
gunzip *.fastq.gz
```
### Download the Annotation File 
  ```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
  ```

#### Gunzip the `Homo_sapiens.GRCh38.113.gtf.gz` file
  ```bash
gunzip Homo_sapiens.GRCh38.113.gtf.gz
  ```

### Move the 'Homo_sapiens.GRCh38.113.gtf' file into the 'annotation.gtf' file
 ```bash
mv Homo_sapiens.GRCh38.113.gtf annotation.gtf
  ```

# Part 3: STAR Genome Indexing

### Step 1: Create STAR Indexing Script
  #### Run:
  ```bash
nano star_indexing.sbatch 
  ```
The `vi` editor can be used here as well.

  #### Paste:
```bash
#!/bin/bash 
#SBATCH --job-name=star_indexing 
#SBATCH --output=star_indexing.out 
#SBATCH --error=star_indexing.err 
#SBATCH --time=4:00:00 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem=100G 
#SBATCH --mail-user=YOUR EMAIL 
#SBATCH --mail-type=ALL 
```
### Submit STAR Indexing Job
```bash
sbatch star_indexing.sbatch
 ```

### You now should have a STAR index
- To check:
  ```bash
  ls -lh star_index
  ```

# Part 3: Read Alignment with STAR

Generate Alignment Script
```bash
nano star_alignment.sbatch
```
Paste, Change Email
```bash
#!/bin/bash
#SBATCH --job-name=STAR_alignment
#SBATCH --output=star_alignment.out
#SBATCH --error=star_alignment.err
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR EMAIL

# Load necessary modules
module load anaconda3/2024.10-1
conda activate bioinfo-env
module load STAR

# Run STAR alignment for each single-end read file
for file in SRR11412215.fastq SRR11412216.fastq SRR11412217.fastq SRR11412218.fastq \
            SRR11412227.fastq SRR11412228.fastq SRR11412229.fastq SRR11412230.fastq SRR11412231.fastq
do
    STAR --genomeDir star_index \
         --readFilesIn $file \
         --outFileNamePrefix star_output/${file%%.fastq}_ \
         --outSAMtype BAM SortedByCoordinate
done
```
### Submit Job
```bash
sbatch star_alignment.sbatch
```
### Verify BAM files exist
```bash
ls -lh star_output/
samtools flagstat star_output/SRR11412215_Aligned.sortedByCoord.out.bam
```
### Part 4: Create Counts File Using Subread

Create Feature Counts Script:
```bash
nano featureCounts.sbatch
```
Paste:
```bash
#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts.out
#SBATCH --error=featureCounts.err
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR EMAIL

# Load necessary modules
module load anaconda3/2024.10-1
conda activate bioinfo-env
module load subread

# Run featureCounts to quantify gene expression
featureCounts -T 8 -a annotation.gtf -o counts.txt star_output/*.bam
```
Submit Job:
```bash
nano featureCounts.sbatch
```
Verify Output
```bash
ls -lh counts.txt
head counts.txt
```
### Part 5 Analysis in R Studio: Setup and Load Required Libraries
#### Ensure that all required packages are installed and loaded:
```r
# Install Bioconductor and DESeq2 if not installed 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 
BiocManager::install("DESeq2") 
 
# Load required libraries 
library(DESeq2) 
library(ggplot2) 
library(pheatmap) 
  ```
### R Studio may need to update these Libraries. Agree by typing in the console and closing and reopening R Studio once updates are complete.
### Load and Prepared Count Data
  - Download Counts File from Student-Led-Tutorial 3
  - Verify File location
  - Set the working directory or use file.choose() to select the count file:
    # Set working directory:
    On Windows:
    # If you do not know the file path use:
```r
print(file.path(Sys.getenv("HOME"), "Downloads"))
```
# Or:
```r
print(file.path(Sys.getenv("USERPROFILE"), "Downloads"))
```
# To set working directory paste in the correct file path:
```r
setwd(file.path(Sys.getenv("USERPROFILE"), "Downloads"))
```
   # For MAC or Linux:
```r
setwd(file.path(Sys.getenv("HOME"), "Downloads"))
```
Example of path:
```r
setwd("C:/Users/19897/Documents/Downloads")
```

# Verify file exists 
```r
file.exists("counts.txt")  # Should return TRUE 
```
### Read and Load count data from the file, again you will have to change the file path, remember to leave the counts.txt file in the path:
  ```r
counts <- read.table("C:/Users/YOUR_USERNAME/Downloads/counts.txt", 
                     header=TRUE, 
                     row.names=1, 
                     sep="\t", 
                     check.names=FALSE)
  ```
### Remove metadata Columns (Chr, Start, End, Strand, Length) 
  ```r
counts <- counts[, -(1:5)]
  ```
### Verify data structure
- `dim(counts)` checks dimensions (genes, samples)
- `colnames(counts)` ensures only sample names remain
```r
dim(counts)   
```
```r
colnames(counts)   
```
# Part 5: Create Sample Metadata
### Define experimental conditions for each sample: 
  ```r
coldata <- data.frame( 
  condition = c(rep("mock", 4), rep("covid", 5)),  # Adjust sample numbers if needed 
  replicate = c(1:4, 1:5)  )  
  ```
### Assign sample names
```r
rownames(coldata) <- colnames(counts)
```
### Check that row names match column names (Response should say: TRUE)
```r
all(rownames(coldata) == colnames(counts))
```
# Part 6: Differential Expression Analysis with DESeq2
### Create DESeq2 dataset
```r
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
```
### Run DESeq2 differential expression analysis
```r
dds <- DESeq(dds)
```
### Get results
```r
results <- results(dds
```
### Save results to CSV
```bash
write.csv(as.data.frame(results), "differential_expression_results.csv")
```
# Part 7: Visualize Results
### Step 1: Volcano Plot
  - Identify Significant genes
```r
results$significant <- !is.na(results$padj) & results$padj < 0.05 & abs(results$log2FoldChange) > 1
```
### Step 2: Generate Volcano Plot
```r
ggplot(results, aes(x=log2FoldChange, y=-log10(padj), color=significant)) + 
  geom_point(alpha=0.7) + 
  theme_minimal() + 
  scale_color_manual(values=c("gray", "red")) + 
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10 Adjusted p-value")
```
### Step 3: Select the top 50 genes based on adjusted p-value
```r
library(pheatmap) 
top_genes <- rownames(results)[order(results$padj, na.last=NA)][1:50]
```
### Step 4: Normalize couns for heatmap visualization
```r
normalized_counts <- counts(dds, normalized=TRUE)
```
### Step 5: Generate heatmap
```r
pheatmap(normalized_counts[top_genes, ], cluster_rows=TRUE, cluster_cols=TRUE, scale="row")
```
