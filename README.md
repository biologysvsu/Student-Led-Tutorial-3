# Student-Led-Tutorial-3
# Task: Tutorial for Gene Expression Analysis Using STAR and Visualization
# Date: March 13th

## **Objective**
Prepare a tutorial for your peers on performing RNA-seq data alignment using **STAR** and analyzing gene expression. The tutorial will cover alignment, quantification, and visualization of gene expression data, resulting in publication-ready visuals.

---

## **Required Software**
1. **STAR**: For RNA-seq alignment.
   - [STAR Manual](https://github.com/alexdobin/STAR)
2. **FeatureCounts**: For gene expression quantification.
   - [FeatureCounts Documentation](http://bioinf.wehi.edu.au/featureCounts/)
3. **R** or **Python**: For data visualization.
   - R packages: `ggplot2`, `pheatmap`.
   - Python libraries: `matplotlib`, `seaborn`.
4. **Genome and Annotation Files**:
   - Reference genome and GTF annotation file from [Ensembl](https://ftp.ensembl.org/pub/).

---

## **Data to Use**
- **Sample Dataset**: Paired-end RNA-seq data from Homo sapiens (e.g., SRR12345678).
  - Download using SRA Toolkit:
    ```bash
    fastq-dump --split-files SRR12345678
    ```
    - This will generate paired-end files:
      - `SRR12345678_1.fastq` and `SRR12345678_2.fastq`.

---

## **Tasks and Deliverables**
### **Part 1: Preparing Reference Files**
1. Download the **reference genome (FASTA)** and **gene annotation file (GTF)**:
   - Example: GRCh38 human genome from [Ensembl](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/).

2. Generate a STAR genome index:
   ```bash
   STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles reference.fasta --sjdbGTFfile annotation.gtf --sjdbOverhang 100

### **Part 2: Aligning Reads**

1. Align paired-end reads to the reference genome using STAR:
   ```bash
   STAR --genomeDir star_index --readFilesIn SRR12345678_1.fastq SRR12345678_2.fastq --outFileNamePrefix aligned_ --outSAMtype BAM SortedByCoordinate

   - Key output files:
     - aligned_Aligned.sortedByCoord.out.bam: Sorted BAM file with alignments.
     - aligned_Log.final.out: STAR alignment statistics.

2. Check alignment statistics:
   - Total mapped reads.
   - Percentage of uniquely mapped reads.

### **Part 3: Gene Expression Quantification**

1. Use FeatureCounts to count reads mapped to genes:
   ```bash
   featureCounts -a annotation.gtf -o counts.txt aligned_Aligned.sortedByCoord.out.bam
  - Key output:
    - *counts.txt*: Gene-level read counts.
2. Normalize the counts (e.g., TPM or FPKM) using R or Python.
   ```R
   library(dplyr)

   # Load the counts data
   counts <- read.table("counts.txt", header = TRUE, row.names = 1)

   # Add gene lengths (in kb) for normalization (from annotation file)
   counts$Length_kb <- counts$Length / 1000

   # Calculate FPKM
   counts$FPKM <- (counts$Counts / counts$Length_kb) / sum(counts$Counts) * 1e6

   # Calculate TPM
   counts$TPM <- counts$FPKM / sum(counts$FPKM) * 1e6

   # Save normalized counts
   write.table(counts, "normalized_counts.txt", sep = "\t", quote = FALSE)


### **Part 4: Data Visualization**

1. Load the normalized counts into R or Python.

2. Create visualizations:
- Heatmap of gene expression:
   ```R
   library(pheatmap)
   data <- read.table("normalized_counts.txt", header=TRUE, row.names=1)
   pheatmap(data, scale="row", clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="complete")
- Volcano Plot for differential expression:
   ```R
   library(ggplot2)
   data <- read.csv("differential_expression.csv")
   ggplot(data, aes(x=log2FoldChange, y=-log10(pvalue), color=significant)) +
      geom_point() +
      theme_minimal()
- Bar Chart:
-   Show expression levels of top 10 differentially expressed genes.

### **Optional: Use IGV to visualize mapped reads**:
- Load the BAM file and gene annotation to check alignment coverage.
