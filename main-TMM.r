# Install from CRAN
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("FactoMineR")
install.packages("reshape2")
install.packages("shiny")
install.packages("ggpubr")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")

# Load necessary libraries
library(limma)
library(edgeR)
library(ggplot2)
library(gridExtra)
library(FactoMineR)
library(reshape2)
library(shiny)
library(ggpubr)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# Load data from csv file
count <- read.csv("countdata.csv", header = T, row.names = 1)

# Convert dataframe to matrix
countmatrix <- as.matrix(count)

# Load sample and batch data from csv file
coldata <- read.csv("coldata.csv", header = T, row.names = 1)
coldata$group <- factor(coldata$group)
coldata$batch <- factor(coldata$batch)

# Create DGEList object
dge <- DGEList(counts = countmatrix)

# Filter the data
keep <- filterByExpr(dge, group = coldata$group)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize counts using TMM method
dge <- calcNormFactors(dge)

# Set "C" group as the reference level ("H"/"C")
coldata$group <- relevel(coldata$group, ref = "C")

# Create design matrix and perform voom transformation
design <- model.matrix(~ batch + group, data = coldata)
v <- voom(dge, design, plot = TRUE)
v_raw <- v$E
write.csv(v_raw, file = "voom_counts.csv")

## DEG analysis (before batch effect removal) ##
# Set an FDR threshold and fold change threshold
fdr_threshold <- 0.05
logFC_threshold <- log2(1.5)

# Fit linear model on raw data
fit_raw <- lmFit(v_raw, design)

# Perform hypothesis testing
fit_raw <- eBayes(fit_raw)

# Extract results
results_raw <- topTable(fit_raw, coef="groupH", number=Inf, sort.by="p")

# Adjust for multiple testing (FDR)
results_raw$adj.P.Val <- p.adjust(results_raw$P.Value, method="BH")

# Select differentially expressed genes
de_genes_raw <- results_raw[(results_raw$adj.P.Val < fdr_threshold) & (abs(results_raw$logFC) > logFC_threshold),]

# Remove batch effect
v_corrected <- removeBatchEffect(v_raw, batch = coldata$batch)

# Save corrected data to a file
write.csv(v_corrected, file = "voom_batch_corrected_counts.csv")

## DEG analysis (after batch effect removal) ##
# Create new design matrix without batch
design_no_batch <- model.matrix(~ group, data = coldata)

# Fit linear model
fit_corrected <- lmFit(v_corrected, design_no_batch)

# Perform hypothesis testing
fit_corrected <- eBayes(fit_corrected)

# Extract results
results_corrected <- topTable(fit_corrected, coef="groupH", number=Inf, sort.by="p")

# Adjust for multiple testing (FDR)
results_corrected$adj.P.Val <- p.adjust(results_corrected$P.Value, method="BH")

# Select differentially expressed genes
de_genes_corrected <- results_corrected[(results_corrected$adj.P.Val < fdr_threshold) & (abs(results_corrected$logFC) > logFC_threshold),]

# Convert ENSG ID to gene symbol for raw data
gene_symbols_raw <- mapIds(org.Hs.eg.db, keys=rownames(de_genes_raw), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
de_genes_raw$gene_symbol <- gene_symbols_raw

# Convert ENSG ID to gene symbol for corrected data
gene_symbols_corrected <- mapIds(org.Hs.eg.db, keys=rownames(de_genes_corrected), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
de_genes_corrected$gene_symbol <- gene_symbols_corrected

# Save the differentially expressed genes with gene symbols to a file
write.csv(de_genes_raw, file="differentially_expressed_genes_raw.csv")
write.csv(de_genes_corrected, file="differentially_expressed_genes_corrected.csv")

# Sort by adjusted P-value and absolute logFC
results_corrected_sorted <- results_corrected[order(results_corrected$adj.P.Val, -abs(results_corrected$logFC)),]
results_raw_sorted <- results_raw[order(results_raw$adj.P.Val, -abs(results_raw$logFC)),]

# Extract the top 100 genes
top100_genes_corrected <- head(results_corrected_sorted, 100)
top100_genes_raw <- head(results_raw_sorted, 100)

# Extract the corresponding expression data for the top genes
top100_gene_data_corrected <- v_corrected[rownames(top100_genes_corrected),]
top100_gene_data_raw <- v_raw[rownames(top100_genes_raw),]

# Save the expression data of top 100 genes to a file
write.csv(top100_gene_data_corrected, file="top100_gene_data_corrected.csv")
write.csv(top100_gene_data_raw, file="top100_gene_data_raw.csv")
## http://www.heatmapper.ca/expression/
# Scale type Row
# Colour Brightness 0
# Numver of Shades 50
# Colour Scheme  #0016DB #FFFFFF #FF1900 Missing Data #000000
# Clustering Method Average Linkage
# Euclidean
# Aplly Clustering Rows
# DPI 300
# Plot Width 600 Height 600

## Volvanoplot ##
# Create volcano plot using EnhancedVolcano
p_corrected <- EnhancedVolcano(results_corrected,
    lab = rownames(results_corrected),
    x = 'logFC',
    y = 'adj.P.Val',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[10]~ 'FDR'),
    xlim = c(-2, 2),
    ylim = c(0, 2), # Set minimum and maximum for y-axis
    title = 'Differentially expressed genes',
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 3.0,
    labSize = 3.0,
    selectLab = character(0)) # Hide all gene labels

# Save pdf and png file
ggsave('volcano_corrected.pdf', plot = p_corrected, width = 10, height = 10, dpi = 300)
ggsave('volcano_corrected.png', plot = p_corrected, width = 10, height = 10, dpi = 300)

p_raw <- EnhancedVolcano(results_raw,
    lab = rownames(results_raw),
    x = 'logFC',
    y = 'adj.P.Val',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[10]~ 'FDR'),
    xlim = c(-2, 2),
    ylim = c(0, 2), # Set minimum and maximum for y-axis
    title = 'Differentially expressed genes',
    pCutoff = 0.05,
    FCcutoff = log2(1.5),
    pointSize = 3.0,
    labSize = 3.0,
    selectLab = character(0)) # Hide all gene labels
# Save pdf and png file
ggsave('volcano_raw.pdf', plot = p_raw, width = 10, height = 10, dpi = 300)
ggsave('volcano_raw.png', plot = p_raw, width = 10, height = 10, dpi = 300)

## Box plot ##
# Transform matrix to data frame and rename columns
v_corrected_df <- as.data.frame(melt(v_corrected))
v_raw_df <- as.data.frame(melt(v_raw))
colnames(v_corrected_df) <- c("Gene", "Sample", "normalized_count")
colnames(v_raw_df) <- c("Gene", "Sample", "normalized_count")

# Create box plot
after_plot <- ggboxplot(v_corrected_df, x = "Sample", y = "normalized_count")
before_plot <- ggboxplot(v_raw_df, x = "Sample", y = "normalized_count")

# Save plot to a pdf file
ggsave("corrected_boxplot.pdf", after_plot, width = 16, height = 9)
ggsave("raw_boxplot.pdf", before_plot, width = 16, height = 9)


## PCA analysis ##
# Perform PCA on the data before batch correction
pca_before_voom <- prcomp(t(v_raw))
pca_before_voom_df <- as.data.frame(pca_before_voom$x[, 1:2])
colnames(pca_before_voom_df) <- c("PC1", "PC2")
pca_before_voom_df$group <- coldata$group
pca_before_voom_df$batch <- coldata$batch

# Plot PCA before batch correction (Voom)
plot_before_voom <- ggplot(pca_before_voom_df, aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA before batch correction (Voom)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "bottom")

# Perform PCA on the data after batch correction
pca_after_voom <- prcomp(t(v_corrected))
pca_after_voom_df <- as.data.frame(pca_after_voom$x[, 1:2])
colnames(pca_after_voom_df) <- c("PC1", "PC2")
pca_after_voom_df$group <- coldata$group
pca_after_voom_df$batch <- coldata$batch

# Plot PCA after batch correction (Voom)
plot_after_voom <- ggplot(pca_after_voom_df, aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA after batch correction (Voom)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "bottom")

write.csv(pca_before_voom_df, file = "pca_before_voom.csv")
write.csv(pca_after_voom_df, file = "pca_after_voom.csv")

ggsave("pca_before_voom.png", plot_before_voom, width = 8, height = 6, units = "in", dpi = 300)
ggsave("pca_after_voom.png", plot_after_voom, width = 8, height = 6, units = "in", dpi = 300)

# Create a grid with both PCA plots
merged_plots <- grid.arrange(ggplotGrob(plot_before_voom), ggplotGrob(plot_after_voom), ncol = 2)

# Save the merged PCA plots to a file
ggsave("merged_PCA_plots.png", plot = merged_plots, width = 16, height = 6)
