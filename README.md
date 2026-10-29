# Student-Led-Tutorial-3
# Task: Tutorial for Gene Expression Analysis Using STAR, Gene Expression Quantification, and Differential Expression Analysis

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
  - SRA Accessions: `SRR11412227`, `SRR11412228`, `SRR11412229`, `SRR11412230` (4 replicates).
There are more data available in this project, students should be free to explore these data and modify SRA accessions accordingly.

## **Tasks and Deliverables**
### **Part 1: Data Preparation**
1. Download all replicates for mock and COVID-19-infected samples:
        Mock-infected: SRR11412215 to SRR11412218.
        COVID-19-infected: SRR11412227 to SRR11412230.
### Step 1. Create a directory for your transcriptome data
0. Request resources:
```
interact -t 4:00:00 --ntasks-per-node=48
```

1. Navigate to your personal folder in the Ocean project directory:

   `cd /ocean/projects/bio260081p/your-username`

   Replace `your-username` with your actual username.

2. Create a new folder named `transcriptome` and move into it:

   `mkdir -p transcriptome`

   `cd transcriptome`

3. Confirm that you are in the correct directory by typing:

   `pwd`

   Your path should look like this:

   `/ocean/projects/bio260081p/your-username/transcriptome`

### Step 2. Create a file containing the SRA accession numbers

1. Open the `vi` text editor to create a file named `srr_accessions.txt`:

   `vi srr_accessions.txt`

2. Press lowercase **`i`** to enter *Insert Mode*.

3. Paste the following accession numbers, **one per line**:

   ```
   SRR11412215
   SRR11412216
   SRR11412217
   SRR11412218
   SRR11412227
   SRR11412228
   SRR11412229
   SRR11412230
   ```

4. Save the file and exit `vi`:
   - Press **Esc** to leave Insert Mode.
   - Type **`:wq`**.
   - Press **Enter** to save and exit.

5. Verify that your file was saved correctly:

   `cat srr_accessions.txt`

   You should see all nine accession numbers displayed, one per line.
```
module load anaconda3
conda create -n sra-tools -c conda-forge -c bioconda sra-tools -y
```

```
vi download_sra.sh
```
- Type I
```

#!/bin/bash
set -e

module load anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sra-tools

mkdir -p sra-data

while read -r SRR; do
    echo "Downloading $SRR..."

    prefetch --max-size 100G "$SRR" -O "$PWD/$SRR"

    fasterq-dump "$PWD/$SRR/$SRR.sra" \
        --threads 8 \
        --split-files \
        -O fastq_files/

    gzip -1 fastq_files/"$SRR"*.fastq

done < srr_accessions.txt

conda deactivate
```
-script may give you an error if anaconda3 was loaded before script was run, exit the HPC and run it again.
-Run script using
```
bash download_sra.sh &
```
3. Download the reference genome (FASTA) and annotation file (GTF):
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
interact -t 4:00:00 --ntasks-per-node=48
```
6. **(OPTIONAL)** Generate the STAR genome index (This step requires 3 hours to complete in 16 cpus with total RAM of 50 Gigabytes):
   ```bash
  mkdir -p star_index 
  
  STAR --runMode genomeGenerate \
  --runThreadN 48 \
  --genomeDir star_index \
  --genomeFastaFiles reference.fasta \
  --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 100

# START HERE

7. Create a symlink to copy index files to your current directory
```
ln -s /ocean/projects/bio260081p/shared/tutorial-data/tutorial_3_data/star_index .
ln -s /ocean/projects/bio260081p/shared/tutorial-data/tutorial_3_data/sra_data .
ln -s /ocean/projects/bio260081p/shared/tutorial-data/tutorial_3_data/reference.fasta .


```
### **Part 2: Align Reads with STAR**
1. Open an interactive session to run this alignment
```
salloc --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=03:00:00 --partition=RM
```
2. Align each sample (mock ) to the genome:
```bash
mkdir mock
cd mock

STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412218.fastq \
     --outFileNamePrefix mock_rep4_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
- Replace SRR11412215 with the appropriate sample name for each replicate (4).
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412216.fastq \
     --outFileNamePrefix mock_rep2_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412217.fastq \
     --outFileNamePrefix mock_rep3_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412218.fastq \
     --outFileNamePrefix mock_rep4_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
2. Repeat with COVID-infected samples:
```bash
cd ..
mkdir covid
cd covid

STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412231.fastq \
     --outFileNamePrefix covid_rep5_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
- Replace SRR11412227 with the appropriate sample name for each replicate (5).
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412228.fastq \
     --outFileNamePrefix covid_rep2_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412229.fastq \
     --outFileNamePrefix covid_rep3_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412230.fastq \
     --outFileNamePrefix covid_rep4_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```
```
STAR --genomeDir ../star_index \
     --readFilesIn ../sra_data/SRR11412231.fastq \
     --outFileNamePrefix covid_rep5_ \
     --runThreadN 8 \
     --outSAMtype BAM SortedByCoordinate
```

3. Output files for each sample:
- Sorted BAM file (e.g., `mock_rep1_Aligned.sortedByCoord.out.bam`).
- Alignment log file (e.g., `mock_rep1_Log.final.out`).

### **Part 3: Gene Expression Quantification**
1.Install `subread`
```
module load anaconda3
```
```
conda create -n bioinfo -c bioconda subread
```
```
conda activate bioinfo
```

2. Use FeatureCounts to count reads mapped to genes for all samples and replicates :
```
cd ..
```
```bash
   featureCounts -T 8 -a annotation.gtf -o counts.txt \
   mock/mock_rep1_Aligned.sortedByCoord.out.bam \
   mock/mock_rep2_Aligned.sortedByCoord.out.bam \
   mock/mock_rep3_Aligned.sortedByCoord.out.bam \
   mock/mock_rep4_Aligned.sortedByCoord.out.bam \
   covid/covid_rep1_Aligned.sortedByCoord.out.bam \
   covid/covid_rep2_Aligned.sortedByCoord.out.bam \
   covid/covid_rep3_Aligned.sortedByCoord.out.bam \
   covid/covid_rep4_Aligned.sortedByCoord.out.bam \
   covid/covid_rep5_Aligned.sortedByCoord.out.bam
```
- Include all BAM files from mock and COVID-infected replicates.

2. Output file:
- `counts.txt`: Contains gene-level counts for all samples.
- Explore the `counts.txt.summary` file
```
cat counts.txt.summary
``` 
### **Part 4: Differential Expression Analysis**
1. Install R and RStudio
2. Open RStudio
3. Create and save new project
4. Open new script.
5. Copy, paste and run this script in R.
```R
###STUDENT LED TUTORIAL 3###
###NAME: ______
###DATE: ______
###CONTACT: _________

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Load data
counts <- read.table("counts.txt", header=TRUE, row.names=1, sep="\t")
counts <- counts[, -(1:5)]  # Now only sample columns remain

coldata <- data.frame(
  condition = factor(c(rep("mock", 4), rep("covid", 5))),  # Convert to factor
  replicate = c(1:4, 1:5)
)

rownames(coldata) <- colnames(counts)

# Check if dimensions match
if (ncol(counts) != nrow(coldata)) {
  stop("Error: counts and coldata dimensions do not match!")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get and save results
results <- results(dds)

# Remove NA values
results <- na.omit(results)

# Save results to CSV
write.csv(results, "differential_expression_results.csv")

# View results structure
head(results)

#Create Volcano Plot
library(ggplot2)

# Define thresholds for significance
padj_threshold <- 0.05  # Adjusted p-value cutoff
lfc_threshold <- 1       # Log2 fold change cutoff

# Add significance categories
results$significance <- "Not Significant"
results$significance[results$padj < padj_threshold & results$log2FoldChange > lfc_threshold] <- "Upregulated"
results$significance[results$padj < padj_threshold & results$log2FoldChange < -lfc_threshold] <- "Downregulated"

# Create the volcano plot
ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")

ggsave("volcano_plot.png", width = 8, height = 6, dpi = 300)

#Check the Number of Upregulated vs. Downregulated Genes
table(results$significance)

#Add gene accessions to graph
install.packages("ggrepel")
library(ggrepel)

# Filter significant genes
top_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3) +  # Add gene labels
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")
#Call top differentially expressed genes
top_genes
```
