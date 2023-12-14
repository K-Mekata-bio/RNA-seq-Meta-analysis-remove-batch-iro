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

BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("sva")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("sva")

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(FactoMineR)
library(reshape2)
library(shiny)
library(ggpubr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(sva)
library(limma)
library(edgeR)

# Load data from csv file
count <- read.csv("countdata.csv", header = T, row.names = 1)

# Convert dataframe to matrix
countmatrix <- as.matrix(count)

# Load sample and batch data from csv file
coldata <- read.csv("coldata.csv", header = T, row.names = 1)
coldata$group <- factor(coldata$group)
coldata$paper <- factor(coldata$paper)

# Create a DGEList object from the count data
dge <- DGEList(counts = countmatrix)

# Apply the filter
keep <- filterByExpr(dge, group = coldata$group)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Convert the filtered DGEList object back to a matrix
filtered_countmatrix <- dge$counts

# Normalize counts using log transformation
logcounts <- log2(filtered_countmatrix + 1)
write.csv(logcounts, file = "logcounts.csv")

## Mean variance trend plot
# Perform voom transformation
v <- voom(filtered_countmatrix, plot=TRUE)

## Remove Batch effect
# Set "C" group as the reference level ("H"/"C")
coldata$group <- relevel(coldata$group, ref = "C")

# Create design matrix and perform voom transformation
design <- model.matrix(~group, data = coldata)

# Remove batch effects with ComBat
adjusted_logcounts <- ComBat(logcounts, batch = coldata$batch, mod = design, par.prior = TRUE, prior.plots = FALSE)

write.csv(adjusted_logcounts, file = "adjusted_logcounts.csv")

## DEG analysis (before batch effect removal) ##
# Set an FDR threshold and fold change threshold
fdr_threshold <- 0.05
logFC_threshold <- log2(1.5) # 1.5-fold change

# Fit linear model
fit_raw <- lmFit(logcounts, design)

# Perform hypothesis testing
fit_raw <- eBayes(fit_raw)

# Extract results CHANGE "coef" name!!!!!!
results_raw <- topTable(fit_raw, coef="groupH", number=Inf, sort.by="p")

# Adjust for multiple testing (FDR)
results_raw$adj.P.Val <- p.adjust(results_raw$P.Value, method="BH")

# Select differentially expressed genes
de_genes_raw <- results_raw[(results_raw$adj.P.Val < fdr_threshold) & (abs(results_raw$logFC) > logFC_threshold),]

## DEG analysis (before batch effect removal) ##
# Fit linear model
fit_corrected <- lmFit(adjusted_logcounts, design)

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

# Convert ENSG ID to gene symbol for raw data
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
top100_gene_data_corrected <- adjusted_logcounts[rownames(top100_genes_corrected),]
top100_gene_data_raw <- logcounts[rownames(top100_genes_raw),]

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

## Volcanoplot ##
# Create volcano plot using EnhancedVolcano
p_corrected <- EnhancedVolcano(results_corrected,
    lab = rownames(results_corrected),
    x = 'logFC',
    y = 'adj.P.Val',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~Log[10]~ 'FDR'),
    xlim = c(-2, 2),
    ylim = c(0, 4),
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
    ylim = c(0, 6),
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
adjusted_logcounts_df <- as.data.frame(reshape2::melt(adjusted_logcounts))
logcounts_df <- as.data.frame(reshape2::melt(logcounts))
colnames(adjusted_logcounts_df) <- c("Gene", "Sample", "normalized_count")
colnames(logcounts_df) <- c("Gene", "Sample", "normalized_count")

# Create box plot
after_plot <- ggpubr::ggboxplot(adjusted_logcounts_df, x = "Sample", y = "normalized_count")
before_plot <- ggpubr::ggboxplot(logcounts_df, x = "Sample", y = "normalized_count")

# Save plot to a pdf file
ggplot2::ggsave("corrected_boxplot.pdf", after_plot, width = 16, height = 9)
ggplot2::ggsave("raw_boxplot.pdf", before_plot, width = 16, height = 9)

## PCA analysis ##
# Perform PCA on the data before batch correction
pca_before_combat <- prcomp(t(logcounts))
pca_before_combat_df <- as.data.frame(pca_before_combat$x[, 1:2])
colnames(pca_before_combat_df) <- c("PC1", "PC2")
pca_before_combat_df$group <- coldata$group
pca_before_combat_df$paper <- coldata$paper

# Plot PCA before batch correction (ComBat)
plot_before_combat <- ggplot(pca_before_combat_df, aes(x = PC1, y = PC2, color = group, shape = paper)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA before batch correction (ComBat)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "bottom")

# Perform PCA on the data after batch correction
pca_after_combat <- prcomp(t(adjusted_logcounts))
pca_after_combat_df <- as.data.frame(pca_after_combat$x[, 1:2])
colnames(pca_after_combat_df) <- c("PC1", "PC2")
pca_after_combat_df$group <- coldata$group
pca_after_combat_df$paper <- coldata$paper

# Plot PCA after batch correction (ComBat)
plot_after_combat <- ggplot(pca_after_combat_df, aes(x = PC1, y = PC2, color = group, shape = paper)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA after batch correction (ComBat)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "bottom")

write.csv(pca_before_combat_df, file = "pca_before_combat.csv")
write.csv(pca_after_combat_df, file = "pca_after_combat.csv")

ggsave("pca_before_combat.png", plot_before_combat, width = 8, height = 6, units = "in", dpi = 300)
ggsave("pca_after_combat.png", plot_after_combat, width = 8, height = 6, units = "in", dpi = 300)

# Create a grid with both PCA plots
merged_plots <- grid.arrange(ggplotGrob(plot_before_combat), ggplotGrob(plot_after_combat), ncol = 2)

# Save the merged PCA plots to a file
ggsave("merged_PCA_plots.png", plot = merged_plots, width = 16, height = 6)
