# ============================================================
# 05_figures.R
# Purpose: Volcano plots, heatmap, and PCA for DEG results
# ============================================================

library(EnhancedVolcano)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(stringr)
library(circlize)  # for colorRamp2 in heatmap coloring

# ============================================================
# LOAD DATA
# ============================================================

# Load DEG results
ec  <- read.csv("../data/processed/deg_Entorhinal_Cortex.csv")
sfg <- read.csv("../data/processed/deg_Superior_Frontal_Gyrus.csv")

# Load normalized expression matrix and clean metadata
norm_expr <- readRDS("../data/processed/norm_expr_matrix.rds")
pdata     <- readRDS("../data/processed/metadata.rds")

# Clean metadata — same steps as 03_limma_de.R
pdata$disease <- pdata$`Disease State:ch1`
pdata$region  <- pdata$`Organ Region:ch1`
pdata$region  <- gsub("hippocampus", "Hippocampus", pdata$region)
pdata$region  <- gsub("Posterior Singulate", "Posterior Cingulate", pdata$region)
pdata$disease <- gsub("Alzheimer's Disease", "AD", pdata$disease)

complete_idx <- !is.na(pdata$disease)            &
  !is.na(pdata$region)             &
  nchar(trimws(pdata$disease)) > 0 &
  nchar(trimws(pdata$region))  > 0

pdata_clean <- pdata[complete_idx, ]
pdata_clean$disease <- str_trim(iconv(pdata_clean$disease,
                                      from = "UTF-8", to = "ASCII//TRANSLIT"))
pdata_clean$region  <- str_trim(iconv(pdata_clean$region,
                                      from = "UTF-8", to = "ASCII//TRANSLIT"))

norm_expr_clean <- norm_expr[, rownames(pdata_clean)]

cat("Data loaded and cleaned ✅\n")

# ============================================================
# PLOT 1 — VOLCANO PLOTS
# ============================================================

# EnhancedVolcano automatically colors and labels genes
# FCcutoff = log2 fold change threshold
# pCutoff  = adjusted p-value threshold
# lab      = gene names to label (rownames of results)

make_volcano <- function(deg_results, title, filename) {
  
  png(filename, width = 1000, height = 900)
  
  print(
    EnhancedVolcano(deg_results,
                    lab            = deg_results$gene,
                    x              = "logFC",
                    y              = "adj.P.Val",
                    title          = title,
                    subtitle       = "AD vs normal | adj.P < 0.05 | |logFC| > 1",
                    pCutoff        = 0.05,
                    FCcutoff       = 1,
                    pointSize      = 2,
                    labSize        = 3.5,
                    col            = c("grey30", "steelblue", "orange", "red"),
                    colAlpha       = 0.6,
                    legendPosition = "right",
                    drawConnectors = TRUE,
                    max.overlaps   = 20)
  )
  
  dev.off()
  cat("Volcano plot saved:", filename, "✅\n")
}

make_volcano(ec,
             "Entorhinal Cortex — AD vs normal",
             "../figures/r_output/11_volcano_EC.png")

make_volcano(sfg,
             "Superior Frontal Gyrus — AD vs normal",
             "../figures/r_output/12_volcano_SFG.png")

# ============================================================
# PLOT 2 — COMPLEXHEATMAP
# ============================================================

# Take top 50 DEGs from EC by absolute logFC
# These are the most strongly differentially expressed genes
top50 <- ec[order(abs(ec$logFC), decreasing = TRUE), ][1:50, ]
top50_genes <- top50$gene

cat("\nTop 50 EC DEGs selected for heatmap\n")

# Subset expression matrix to top 50 genes and all 105 samples
heatmap_mat <- as.matrix(norm_expr_clean[top50_genes, ])

# Scale each gene (row) to mean=0, sd=1
# This makes patterns visible — without scaling, highly expressed
# genes dominate and low expressed genes are invisible
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# --- Build annotation for columns (samples) ---
# This adds colored bars on top of the heatmap showing
# which samples are AD vs normal and which brain region

# Color mappings
disease_colors <- c("AD"     = "#E74C3C",   # red for AD
                    "normal" = "#3498DB")    # blue for normal

region_colors <- c("Entorhinal Cortex"     = "#E8B4B8",
                   "Hippocampus"           = "#A8D5A2",
                   "Medial Temporal Gyrus" = "#85C1E9",
                   "Superior Frontal Gyrus"= "#F9E79F",
                   "Primary Visual Cortex" = "#D7BDE2",
                   "Posterior Cingulate"   = "#FDEBD0")

# HeatmapAnnotation creates the colored bars at the top
col_annotation <- HeatmapAnnotation(
  Disease = pdata_clean$disease,
  Region  = pdata_clean$region,
  col     = list(Disease = disease_colors,
                 Region  = region_colors),
  annotation_name_side = "left"
)

# --- Color scale for expression values ---
# colorRamp2 maps expression values to colors
# Negative = blue (low), Zero = white, Positive = red (high)
col_fun <- colorRamp2(c(-2, 0, 2),
                      c("steelblue", "white", "firebrick"))

# --- Draw heatmap ---
png("../figures/r_output/13_heatmap_top50.png",
    width = 1200, height = 900)

draw(
  Heatmap(heatmap_mat_scaled,
          name                  = "Z-score",
          col                   = col_fun,
          top_annotation        = col_annotation,
          show_column_names     = FALSE,  # too many samples to label
          show_row_names        = TRUE,
          row_names_gp          = gpar(fontsize = 8),
          cluster_rows          = TRUE,   # cluster genes
          cluster_columns       = TRUE,   # cluster samples
          column_title          = "Top 50 DEGs — Entorhinal Cortex (AD vs normal)",
          column_title_gp       = gpar(fontsize = 14, fontface = "bold"),
          heatmap_legend_param  = list(title = "Z-score"))
)

dev.off()
cat("Heatmap saved ✅\n")

# Save the top50 matrix as CSV for Python and JS later
write.csv(as.data.frame(heatmap_mat_scaled),
          "../data/processed/top50_matrix.csv")
cat("top50_matrix.csv saved ✅\n")

# ============================================================
# PLOT 3 — UPDATED PCA (clean 105 samples)
# ============================================================

# Run PCA on the clean 105 sample matrix
# t() transposes so samples are rows, genes are columns
pca <- prcomp(t(norm_expr_clean), scale. = TRUE)

pca_df <- data.frame(
  PC1     = pca$x[, 1],
  PC2     = pca$x[, 2],
  disease = pdata_clean$disease,
  region  = pdata_clean$region
)

var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

png("../figures/r_output/14_pca_clean.png",
    width = 900, height = 700)

print(
  ggplot2::ggplot(pca_df, ggplot2::aes(x      = PC1,
                                       y      = PC2,
                                       color  = disease,
                                       shape  = region)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("AD"     = "#E74C3C",
                                           "normal" = "#3498DB")) +
    ggplot2::labs(title  = "PCA — GSE5281 clean samples (n=105)",
                  x      = paste0("PC1 (", var_explained[1], "%)"),
                  y      = paste0("PC2 (", var_explained[2], "%)"),
                  color  = "Disease Status",
                  shape  = "Brain Region") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")
)

dev.off()
cat("PCA plot saved ✅\n")

cat("\n=== 05_figures.R complete! ===\n")
cat("Check figures/r_output/ for all 14 PNG files\n")