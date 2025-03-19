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
## **Software and Manuals**
### **Required Software**
1. **STAR**: Spliced Transcript Alignment to a Reference.
   - [STAR GitHub Repository](https://github.com/alexdobin/STAR)
   - [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
2. **SRA Toolkit**: To download RNA-seq data from SRA.
   - [SRA Toolkit Documentation](https://github.com/ncbi/sra-tools)
3. **FastQC**: For quality control of RNA-seq reads.
   - [FastQC Website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. **Samtools**: For sequence data manipulation.
   - [Samtools Documentation](http://www.htslib.org/doc/)
5. **Subread/FeatureCounts**: For gene-level quantification.
   - [Subread User Manual](http://bioinf.wehi.edu.au/subread-package/)
6. **R and RStudio**: For differential expression analysis and visualization.
   - [R Project](https://www.r-project.org/)
   - [RStudio Website](https://posit.co/downloads/)

---
## **Dataset Description**
The dataset comes from the **NCBI BioProject PRJNA615032**, which includes RNA-seq data for human cell lines:
- **Mock-infected cells**:
  - SRA Accessions: `SRR11412215`, `SRR11412216`, `SRR11412217`, `SRR11412218` (4 replicates).
- **COVID-19-infected cells**:
  - SRA Accessions: `SRR11412227`, `SRR11412228`, `SRR11412229`, `SRR11412230`, `SRR11412231` (5 replicates).
There are more data available in this project, students should be free to explore these data and modify SRA accessions accordingly.

### **Download Instructions**
Use the SRA Toolkit (must be installed beforehand if run locally, otherwise available in most HPCs) to download paired-end FASTQ files:

   ``` bash
   # Example for a mock-infected sample. More replicates are always better, so repeat step for each SRA    accession.

prefetch --max-size 100G SRR11412215 #Handles large files efficiently (downloads in chunks to avoid corruption) 
fastq-dump --gzip SRR11412215
```
or

   ``` bash
   # Example for a COVID-19-infected sample. More replicates are always better, so repeat step for each SRA accession.
prefetch --max-size 100G SRR11412227
fastq-dump --gzip SRR11412227
```

## **Tasks and Deliverables**
### **Part 1: Data Preparation**
1. Download all replicates for mock and COVID-19-infected samples:
        Mock-infected: SRR11412215 to SRR11412218.
        COVID-19-infected: SRR11412227 to SRR11412231.
2. Download the reference genome (FASTA) and annotation file (GTF):
Source: Ensembl GRCh38 human genome:
- FASTA: Download FASTA file:
```
# Download the fasta file
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename it to reference.fasta
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta
```

- GTF: Download GTF file:
```
# Download the GTF file
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

# Unzip the GTF file
gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Rename it to annotation.gtf
mv Homo_sapiens.GRCh38.113.gtf annotation.gtf
```

3. Load STAR on the HPC
  ```
  module load STAR
  ```
4. **(OPTIONAL)** Run an interactive session requesting high memory resources
```
salloc --mem=50G --time=4:00:00 --cpus-per-task=16
```
6. **(OPTIONAL)** Generate the STAR genome index (This step requires 3 hours to complete in 16 cpus with total RAM of 50 Gigabytes):
   ```bash
   STAR --runMode genomeGenerate --genomeDir star_index \
    --genomeFastaFiles reference.fasta --sjdbGTFfile annotation.gtf --sjdbOverhang 100
7. Create a symlink to copy index files to your current directory
```
ln -s /ocean/projects/agr250001p/shared/tutorial-data/tutorial_7_data/star_index .
```
### **Part 2: Align Reads with STAR**
1. Open an interactive session to run this alignment

2. Align each sample (mock or COVID-infected) to the genome:
   ```bash
STAR --genomeDir star_index \
     --readFilesIn SRR11412215.fastq \
     --outFileNamePrefix mock_rep1_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
- Replace SRR11412215 with the appropriate sample name for each replicate.
3. Output files for each sample:
- Sorted BAM file (e.g., `mock_rep1_Aligned.sortedByCoord.out.bam`).
- Alignment log file (e.g., `mock_rep1_Log.final.out`).

### **Part 3: Gene Expression Quantification**
1. Use FeatureCounts to count reads mapped to genes for all samples and replicates :
   ```bash
   featureCounts -a annotation.gtf -o counts.txt \
   mock_rep1_Aligned.sortedByCoord.out.bam \
   mock_rep2_Aligned.sortedByCoord.out.bam \
   covid_rep1_Aligned.sortedByCoord.out.bam \
   covid_rep2_Aligned.sortedByCoord.out.bam

- Include all BAM files from mock and COVID-infected replicates.

2. Output file:
- `counts.txt`: Contains gene-level counts for all samples.

### **Part 4: Differential Expression Analysis**
1. Use DESeq2 in R for differential expression analysis:
- Load counts.txt and assign sample metadata:
   ```R
   Library(DESeq2)

    # Load data
    counts <- read.table("counts.txt", header=TRUE, row.names=1)
    coldata <- data.frame(
      condition = c(rep("mock", 2), rep("covid", 2)), # Adjust # of reps if needed
      replicate = c(1:2, 1:2)
    )
    rownames(coldata) <- colnames(counts)[-c(1:6)]  # Adjust if needed

    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

    # Perform differential expression
    dds <- DESeq(dds)
    results <- results(dds)
    write.csv(results, "differential_expression_results.csv")
2. Identify differentially expressed genes:

- Filter by adjusted p-value and log2 fold-change (e.g., p.adj < 0.05 and |log2FC| > 1).

### **Part 5: Data Visualization**

1. Volcano Plot:
- Highlight significant genes:
   ```R
   library(ggplot2)
results$significant <- results$padj < 0.05 & abs(results$log2FoldChange) > 1
ggplot(results, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point() + theme_minimal() +
  scale_color_manual(values=c("gray", "red"))
2. Heatmap of Gene Expression:
- Generate a heatmap for top 50 differentially expressed genes:
   ```R
   library(pheatmap)
top_genes <- head(order(results$padj), 50)
normalized_counts <- counts(dds, normalized=TRUE)
pheatmap(normalized_counts[top_genes, ])

### **Optional: Use IGV to visualize mapped reads**:
- Load the BAM file and gene annotation to check alignment coverage.
