## Detailed Summary of the Bioinformatics Data Analysis Script

### Notes
[日本語版](https://github.com/K-Mekata-bio/RNA-seq-Meta-analysis-remove-batch-iro/blob/main/jpabst.md)
- **This code is Version 1. Please modify as needed:**
  - I might make it more automation-friendly if I feel like it. Too hard? Good luck (`･ω･´)
  - I have prepared two types of files, one normalized using the TMM method and the other using log transformation. Please feel free to use whichever you prefer.
- **When removing batch effects, decide which group to use as the reference:**
  - Set "C" group as the reference level ("H"/"C")  coldata$group <- relevel(coldata$group, ref = "C")
- **Please change the threshold values used for DEG identification:**
  - Defaults are, fdr_threshold <- 0.05  logFC_threshold <- log2(1.5).
- **Change the name in coef="" as per your results:**
  - Default is, results_raw <- topTable(fit_raw, coef="groupH", number=Inf, sort.by="p").
- **This page was summarized by ChatGPT-4 and edited:**

### Installation and Library Loading
- **Package Installation from CRAN:**
  - Install `ggplot2`, `gridExtra`, `FactoMineR`, `reshape2`, `shiny`, `ggpubr`
- **Package Installation from Bioconductor:**
  - Install `EnhancedVolcano`, `org.Hs.eg.db`, `sva`, `edgeR`, `limma`
- **Loading Libraries:**
  - Load all the packages installed previously

### Data Loading and Processing
- **Loading Gene Expression Data (`countdata.csv`) and Sample Metadata (`coldata.csv`):**
  - Load from CSV files with headers, setting the first column as row names
- **Data Preprocessing:**
  - Convert gene expression data to a matrix
  - Create and filter DGE list objects
  - Normalize with Log transformation and save to CSV

### Removing Batch Effects
- **Using `ComBat` function for batch effect removal:**
  - Remove batch effects from log-transformed data
  - Save log-counts to CSV

### Differential Expression Gene (DEG) Analysis
- **DEG Analysis before and after batch effect removal:**
  - Set FDR and fold change thresholds
  - Apply linear model, hypothesis testing, and extract results
  - Calculate adjusted P-values after multiple testing and select DEGs
  - Map to gene symbols and save to CSV

### Visualization
- **Volcano Plot (using `EnhancedVolcano`):**
  - Create plots with logfoldchange and FDR as axes
  - Save in PDF and PNG formats
- **Box Plot (using `ggpubr`):**
  - Plot normalized counts for each sample
  - Save in PDF format
- **PCA Analysis (using `prcomp`):**
  - Apply PCA analysis to data before batch effect removal

### Data Saving
- Save analysis and visualization results in CSV, PDF, PNG files
