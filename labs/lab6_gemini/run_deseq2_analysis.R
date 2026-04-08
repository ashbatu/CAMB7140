# Load necessary libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(rtracklayer) # For importing GTF
library(GenomicFeatures) # For making TxDb object

# --- Configuration ---
# Path to kallisto output directories
kallisto_output_path <- "kallisto_output"

# Path to GTF annotation file
gtf_file <- "references/gencode.v49.annotation.gtf.gz"

# Sample information (matching studyDesign.txt)
samples <- data.frame(
  sample = c("HS01", "HS02", "HS03", "CL08", "CL10", "CL11"),
  sra_accession = c("SRR8668755", "SRR8668756", "SRR8668757", "SRR8668769", "SRR8668771", "SRR8668772"),
  group = c("healthy", "healthy", "healthy", "disease", "disease", "disease")
)

# Ensure 'group' is a factor with levels for DESeq2
samples$group <- factor(samples$group, levels = c("healthy", "disease"))

# Paths to abundance.tsv files
files <- file.path(kallisto_output_path, samples$sra_accession, "abundance.tsv")
names(files) <- samples$sample

# Check if all files exist
if (!all(file.exists(files))) {
  stop("One or more kallisto abundance.tsv files not found. Please check paths.")
}

# Create tx2gene mapping from GTF
# This can take a few minutes for large GTF files
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# --- Import kallisto data with tximport ---
txi <- tximport(files, type = "kallisto", txOut = FALSE, countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# --- Create DESeq2 object ---
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ group)

# Pre-filter low count genes (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# --- Differential Expression Analysis ---
res <- results(dds, contrast = c("group", "disease", "healthy"))
res_ordered <- res[order(res$padj),]

# Output full results table
write.csv(as.data.frame(res_ordered), file = "tables/deseq2_full_results.csv")

# Output significant results table (e.g., padj < 0.05)
res_sig <- subset(res_ordered, padj < 0.05)
write.csv(as.data.frame(res_sig), file = "tables/deseq2_significant_results.csv")

# --- Plotting ---

# PCA Plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  theme_bw()

ggsave("plots/pca_plot.pdf", pca_plot, width = 7, height = 6)
ggsave("plots/pca_plot.jpg", pca_plot, width = 7, height = 6)
ggsave("plots/pca_plot.png", pca_plot, width = 7, height = 6)

# Volcano Plot
# Add a column for significance and direction
res_df <- as.data.frame(res) %>%
  mutate(
    significance = case_when(
      padj < 0.01 & log2FoldChange > 0 ~ "Significant (Up)",
      padj < 0.01 & log2FoldChange < 0 ~ "Significant (Down)",
      TRUE ~ "Not significant"
    ),
    gene_label = ifelse(padj < 0.01 & abs(log2FoldChange) > 1, rownames(.), NA) # Label top genes
  )

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Significant (Up)" = "red", "Significant (Down)" = "blue", "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = gene_label), na.rm = TRUE, max.overlaps = 50, size = 3) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plots/volcano_plot.pdf", volcano_plot, width = 8, height = 7)
ggsave("plots/volcano_plot.jpg", volcano_plot, width = 8, height = 7)
ggsave("plots/volcano_plot.png", volcano_plot, width = 8, height = 7)

# Heatmap of top 50 differentially expressed genes
top_genes <- head(res_ordered, 50)
mat <- assay(vsd)[rownames(top_genes), ]
mat <- mat - rowMeans(mat) # Center rows

# Annotate columns with group info
df <- as.data.frame(colData(vsd)[, "group", drop = FALSE])

heatmap_plot <- pheatmap(mat, annotation_col = df, show_rownames = FALSE, 
                         main = "Heatmap of Top 50 Differentially Expressed Genes")

# To save pheatmap plots, you typically need to open a device first
pdf("plots/heatmap_plot.pdf", width = 8, height = 10)
print(heatmap_plot)
dev.off()

jpeg("plots/heatmap_plot.jpg", width = 800, height = 1000)
print(heatmap_plot)
dev.off()

png("plots/heatmap_plot.png", width = 800, height = 1000)
print(heatmap_plot)
dev.off()

message("Differential expression analysis and plotting complete.")
