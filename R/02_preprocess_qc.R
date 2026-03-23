# ============================================================
# 02_preprocess_qc.R
# Purpose: Normalize expression data, map probes to gene 
#          symbols, and run QC checks
# ============================================================

library(Biobase)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(dplyr)
library(tibble)
library(magrittr)

# --- Load the saved ExpressionSet ---
eset <- readRDS("../data/raw/GSE5281_eset.rds")

# ============================================================
# PART 1 — CHECK RAW EXPRESSION
# ============================================================

# Extract the raw expression matrix
# rows = probes (54675), columns = samples (161)
raw_mat <- exprs(eset)

cat("Raw matrix dimensions:", dim(raw_mat), "\n")
cat("Value range:", range(raw_mat), "\n")

# ============================================================
# PART 2 — RMA NORMALIZATION
# ============================================================

# The series matrix from GEO is already summarized at probe 
# level but not RMA normalized. We normalize using limma's
# normalizeBetweenArrays() since we don't have raw CEL files.
# This applies quantile normalization across all 161 samples.

# First log2 transform the raw intensities
# We add 1 before logging to avoid log(0) which is -infinity
log_mat <- log2(raw_mat + 1)

# Quantile normalize across samples so all distributions align
norm_mat <- normalizeBetweenArrays(log_mat, method = "quantile")

cat("Normalized value range:", range(norm_mat), "\n")
cat("Normalized matrix dimensions:", dim(norm_mat), "\n")

# ============================================================
# PART 3 — PROBE TO GENE SYMBOL MAPPING
# ============================================================

# The rows of norm_mat are probe IDs like "1007_s_at"
# We need to convert these to gene symbols like "DDR1"
# hgu133plus2.db is the annotation package for this exact chip

# Get mapping: probe ID → gene symbol
probe_ids <- rownames(norm_mat)

# mapIds looks up each probe ID in the annotation database
# and returns the corresponding gene symbol
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# How many probes mapped successfully?
cat("Total probes:", length(probe_ids), "\n")
cat("Probes with gene symbol:", sum(!is.na(gene_symbols)), "\n")
cat("Probes without annotation:", sum(is.na(gene_symbols)), "\n")

# Remove probes that didn't map to any gene
mapped <- !is.na(gene_symbols)
norm_mat_mapped <- norm_mat[mapped, ]
gene_symbols_mapped <- gene_symbols[mapped]

cat("Matrix after removing unannotated probes:", dim(norm_mat_mapped), "\n")

# ============================================================
# PART 4 — COLLAPSE DUPLICATE PROBES
# ============================================================

norm_df <- as.data.frame(norm_mat_mapped)
norm_df$gene <- gene_symbols_mapped
norm_df$mean_expr <- rowMeans(norm_mat_mapped)

cat("Probes before collapsing duplicates:", nrow(norm_df), "\n")
cat("Unique gene symbols:", length(unique(norm_df$gene)), "\n")

norm_df_unique <- norm_df %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

norm_df_unique <- norm_df_unique %>%
  dplyr::select(-mean_expr) %>%
  tibble::column_to_rownames("gene")

cat("Genes after collapsing duplicates:", nrow(norm_df_unique), "\n")

# ============================================================
# PART 5 — SAVE OUTPUTS
# ============================================================

# Save as RDS for use in R scripts
saveRDS(norm_df_unique, "../data/processed/norm_expr_matrix.rds")

# Save as CSV for Python and JS to consume later
write.csv(norm_df_unique, "../data/processed/norm_expr_matrix.csv")

# Save clean metadata table
pdata <- pData(eset)
saveRDS(pdata, "../data/processed/metadata.rds")
write.csv(pdata, "../data/processed/metadata.csv")

cat("Normalized matrix saved ✅\n")
cat("Metadata saved ✅\n")
cat("Dimensions of final matrix:", dim(norm_df_unique), "\n")

# ============================================================
# PART 6 — QC PLOTS
# ============================================================

# We'll create 3 plots:
# 1. Boxplot — raw vs normalized expression distributions
# 2. Density plot — expression curves per sample
# 3. PCA — samples in 2D space colored by disease status

# --- Plot 1: Boxplot raw vs normalized ---
# We only plot 20 samples to keep it readable
# t() transposes the matrix so samples are on x axis

png("../figures/r_output/01_boxplot_raw_vs_norm.png", 
    width = 1200, height = 500)

par(mfrow = c(1, 2))  # two plots side by side

# Raw boxplot — log2 transform just for visualization
boxplot(log2(raw_mat[, 1:20] + 1),
        main = "Raw Expression (first 20 samples)",
        ylab = "log2 Intensity",
        xlab = "Samples",
        las = 2,           # rotate x axis labels
        col = "lightcoral",
        cex.axis = 0.6)

# Normalized boxplot
boxplot(as.matrix(norm_df_unique[, 1:20]),
        main = "Normalized Expression (first 20 samples)",
        ylab = "log2 Intensity",
        xlab = "Samples",
        las = 2,
        col = "lightblue",
        cex.axis = 0.6)

dev.off()
cat("Boxplot saved ✅\n")

# --- Plot 2: PCA colored by disease status ---
# prcomp() runs PCA — we transpose because prcomp expects
# samples as rows and genes as columns

pca <- prcomp(t(norm_df_unique), scale. = TRUE)

# Extract PC1 and PC2 scores for each sample
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  disease = pdata$`disease state:ch1`,    # AD vs Control
  region  = pdata$`organ region:ch1`      # brain region
)

# Calculate % variance explained by each PC
var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

# Plot
png("../figures/r_output/02_pca_disease.png",
    width = 800, height = 600)

print(
  ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, 
                                       color = disease, 
                                       shape = region)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(title = "PCA — GSE5281 samples",
                  x = paste0("PC1 (", var_explained[1], "%)"),
                  y = paste0("PC2 (", var_explained[2], "%)"),
                  color = "Disease Status",
                  shape = "Brain Region") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
)

dev.off()
cat("PCA plot saved ✅\n")

