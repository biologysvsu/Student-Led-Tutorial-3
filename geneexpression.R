###STUDENT LED TUTORIAL 3###
###NAME: ______
###DATE: ______
###CONTACT: _________

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2", ask = FALSE)
}
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
volcano <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")
#Show volcano plot in plots window
print (volcano)

#Save volcano plot to session directory
ggsave("volcano_plot.png", width = 8, height = 6, dpi = 300)

#Check the Number of Upregulated vs. Downregulated Genes
table(results$significance)

#Add gene accessions to graph
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

# Filter significant genes
top_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)

volcanolabeled <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3) +  # Add gene labels
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme(legend.position = "top")

#show volcano in the plots window
print (volcanolabeled)

#Save labeled volcano plot to session directory
ggsave("volcano_plot_labels.png", width = 8, height = 6, dpi = 300)


#Call top differentially expressed genes
top_genes

#Choose one gene accession and search them on https://useast.ensembl.org/index.html
