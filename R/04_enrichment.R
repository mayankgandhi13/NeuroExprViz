# ============================================================
# 04_enrichment.R
# Purpose: GO ORA and GSEA enrichment analysis on DEG results
#          from Entorhinal Cortex and Superior Frontal Gyrus
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(fgsea)
library(dplyr)
library(stringr)
library(ggridges)

# ============================================================
# STEP 1 — LOAD DEG RESULTS
# ============================================================

# Load the two regions that had significant DEGs
ec  <- read.csv("../data/processed/deg_Entorhinal_Cortex.csv")
sfg <- read.csv("../data/processed/deg_Superior_Frontal_Gyrus.csv")

cat("EC DEGs loaded:", nrow(ec), "genes\n")
cat("SFG DEGs loaded:", nrow(sfg), "genes\n")

# ============================================================
# STEP 2 — GENE SYMBOL TO ENTREZ ID CONVERSION
# ============================================================

# clusterProfiler works with Entrez IDs internally
# mapIds converts gene symbols → Entrez IDs using org.Hs.eg.db
# org.Hs.eg.db is the human gene annotation database

convert_to_entrez <- function(gene_symbols) {
  # mapIds returns a named vector: name = symbol, value = entrez ID
  entrez <- mapIds(org.Hs.eg.db,
                   keys      = gene_symbols,
                   column    = "ENTREZID",
                   keytype   = "SYMBOL",
                   multiVals = "first")
  # Remove NAs — genes that couldn't be mapped
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

# Convert all genes (background) and significant DEGs for each region
# Background = all genes we tested in limma
all_genes_entrez <- convert_to_entrez(ec$gene)

# Significant DEGs per region
ec_sig  <- ec[ec$adj.P.Val < 0.05 & abs(ec$logFC) > 1, ]
sfg_sig <- sfg[sfg$adj.P.Val < 0.05 & abs(sfg$logFC) > 1, ]

cat("\nEC significant DEGs:", nrow(ec_sig), "\n")
cat("SFG significant DEGs:", nrow(sfg_sig), "\n")

# Convert significant DEGs to Entrez IDs
ec_sig_entrez  <- convert_to_entrez(ec_sig$gene)
sfg_sig_entrez <- convert_to_entrez(sfg_sig$gene)

cat("EC Entrez IDs mapped:", length(ec_sig_entrez), "\n")
cat("SFG Entrez IDs mapped:", length(sfg_sig_entrez), "\n")

# ============================================================
# STEP 3 — GO ORA (Over-Representation Analysis)
# ============================================================

# We run ORA separately for each region
# ont = "BP" means Biological Process only
# pAdjustMethod = "BH" is Benjamini-Hochberg FDR correction
# pvalueCutoff = 0.05 keeps only significant terms
# qvalueCutoff = 0.2 is a secondary FDR filter
# universe = all genes we tested in limma (the background)

run_ora <- function(sig_entrez, background_entrez, label) {
  
  cat("\nRunning GO ORA for:", label, "\n")
  
  ora <- enrichGO(gene          = as.character(sig_entrez),
                  universe      = as.character(background_entrez),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)  # converts Entrez IDs back to gene symbols
  
  cat("Enriched GO terms found:", nrow(ora@result[ora@result$p.adjust < 0.05, ]), "\n")
  return(ora)
}

# Run ORA for both regions
ora_ec  <- run_ora(ec_sig_entrez,  all_genes_entrez, "Entorhinal Cortex")
ora_sfg <- run_ora(sfg_sig_entrez, all_genes_entrez, "Superior Frontal Gyrus")

# Save results as CSV
write.csv(ora_ec@result,  "../data/processed/ora_EC.csv",  row.names = FALSE)
write.csv(ora_sfg@result, "../data/processed/ora_SFG.csv", row.names = FALSE)
cat("\nORA results saved ✅\n")

# ============================================================
# STEP 3B — ORA PLOTS
# ============================================================

# --- Dotplot ---
# x axis = gene ratio (fraction of your DEGs in that GO term)
# dot size = number of genes
# dot color = adjusted p-value (darker = more significant)

png("../figures/r_output/03_ora_dotplot_EC.png",
    width = 900, height = 800)
print(dotplot(ora_ec, showCategory = 20,
              title = "GO Biological Process — Entorhinal Cortex"))
dev.off()

png("../figures/r_output/04_ora_dotplot_SFG.png",
    width = 900, height = 800)
print(dotplot(ora_sfg, showCategory = 20,
              title = "GO Biological Process — Superior Frontal Gyrus"))
dev.off()

cat("ORA dotplots saved ✅\n")

# --- Barplot ---
# bar length = gene count in that GO term
# color = adjusted p-value

png("../figures/r_output/05_ora_barplot_EC.png",
    width = 900, height = 700)
print(barplot(ora_ec, showCategory = 15,
              title = "Top GO Terms — Entorhinal Cortex"))
dev.off()

png("../figures/r_output/06_ora_barplot_SFG.png",
    width = 900, height = 700)
print(barplot(ora_sfg, showCategory = 15,
              title = "Top GO Terms — Superior Frontal Gyrus"))
dev.off()

cat("ORA barplots saved ✅\n")

# ============================================================
# STEP 4 — PREPARE RANKED GENE LISTS FOR GSEA
# ============================================================

# GSEA needs a named numeric vector sorted by logFC
# name = gene symbol, value = logFC
# Sorted descending — most upregulated genes at top

prepare_ranked_list <- function(deg_results) {
  # Sort by logFC descending
  deg_sorted <- deg_results[order(deg_results$logFC, decreasing = TRUE), ]
  
  # Create named vector
  ranked <- deg_sorted$logFC
  names(ranked) <- deg_sorted$gene
  
  # Remove duplicates if any
  ranked <- ranked[!duplicated(names(ranked))]
  
  return(ranked)
}

ranked_ec  <- prepare_ranked_list(ec)
ranked_sfg <- prepare_ranked_list(sfg)

cat("EC ranked list length:", length(ranked_ec), "\n")
cat("SFG ranked list length:", length(ranked_sfg), "\n")

# ============================================================
# STEP 5 — GSEA USING HALLMARK GENE SETS
# ============================================================

# We use the Hallmark gene sets from MSigDB
# These are 50 well-defined biological pathways
# Much more interpretable than raw GO terms
# msigdbr package provides these gene sets in R

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}
library(msigdbr)

# Get Hallmark gene sets for Homo sapiens
# Returns a dataframe with gene set name and gene symbol columns
hallmark_sets <- msigdbr(species = "Homo sapiens", 
                         category = "H")

# Convert to the list format fgsea expects
# fgsea needs a named list: name = pathway, value = vector of gene symbols
hallmark_list <- split(hallmark_sets$gene_symbol, 
                       hallmark_sets$gs_name)

cat("Hallmark gene sets loaded:", length(hallmark_list), "\n")

# Run GSEA for both regions
# minSize/maxSize filter out gene sets that are too small or too large
# nPermSimple = number of permutations for p-value estimation

run_gsea <- function(ranked_list, label) {
  cat("\nRunning GSEA for:", label, "\n")
  
  set.seed(42)  # for reproducibility
  
  gsea_result <- fgsea(pathways  = hallmark_list,
                       stats     = ranked_list,
                       minSize   = 15,
                       maxSize   = 500,
                       nPermSimple = 1000)
  
  # Sort by NES (Normalized Enrichment Score)
  # Positive NES = pathway upregulated in AD
  # Negative NES = pathway downregulated in AD
  gsea_result <- gsea_result[order(gsea_result$NES, decreasing = TRUE), ]
  
  cat("Significant pathways (padj < 0.05):", 
      sum(gsea_result$padj < 0.05, na.rm = TRUE), "\n")
  
  return(gsea_result)
}

gsea_ec  <- run_gsea(ranked_ec,  "Entorhinal Cortex")
gsea_sfg <- run_gsea(ranked_sfg, "Superior Frontal Gyrus")

# Save results — convert leadingEdge list column to string first
# leadingEdge is a list column which can't be written to CSV directly
gsea_ec$leadingEdge  <- sapply(gsea_ec$leadingEdge,  paste, collapse = ";")
gsea_sfg$leadingEdge <- sapply(gsea_sfg$leadingEdge, paste, collapse = ";")

write.csv(gsea_ec,  "../data/processed/gsea_EC.csv",  row.names = FALSE)
write.csv(gsea_sfg, "../data/processed/gsea_SFG.csv", row.names = FALSE)
cat("\nGSEA results saved ✅\n")

# ============================================================
# STEP 5B — GSEA PLOTS
# ============================================================

# --- Ridge plot ---
# Shows logFC distribution of genes in each top pathway
# Peaks to the right = pathway genes tend to be upregulated in AD
# Peaks to the left = pathway genes tend to be downregulated in AD

# Re-run GSEA using clusterProfiler's gseGO for ridge/dot plots
# fgsea is faster but enrichplot plots work better with clusterProfiler objects

run_gsea_cp <- function(ranked_list, label) {
  cat("\nRunning clusterProfiler GSEA for:", label, "\n")
  
  # Convert gene symbols to Entrez IDs for clusterProfiler
  entrez_map <- mapIds(org.Hs.eg.db,
                       keys    = names(ranked_list),
                       column  = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Keep only mapped genes
  mapped      <- !is.na(entrez_map)
  ranked_entrez        <- ranked_list[mapped]
  names(ranked_entrez) <- entrez_map[mapped]
  
  # Remove duplicates
  ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
  
  gsea_cp <- gseGO(geneList     = ranked_entrez,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 15,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
  
  return(gsea_cp)
}

gsea_cp_ec  <- run_gsea_cp(ranked_ec,  "Entorhinal Cortex")
gsea_cp_sfg <- run_gsea_cp(ranked_sfg, "Superior Frontal Gyrus")

# Ridge plot
png("../figures/r_output/07_gsea_ridge_EC.png",
    width = 1000, height = 800)
print(ridgeplot(gsea_cp_ec) + 
        ggplot2::labs(title = "GSEA Ridge Plot — Entorhinal Cortex"))
dev.off()

png("../figures/r_output/08_gsea_ridge_SFG.png",
    width = 1000, height = 800)
print(ridgeplot(gsea_cp_sfg) + 
        ggplot2::labs(title = "GSEA Ridge Plot — Superior Frontal Gyrus"))
dev.off()

# Dotplot
png("../figures/r_output/09_gsea_dot_EC.png",
    width = 900, height = 800)
print(dotplot(gsea_cp_ec, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggplot2::labs(title = "GSEA Dotplot — Entorhinal Cortex"))
dev.off()

png("../figures/r_output/10_gsea_dot_SFG.png",
    width = 900, height = 800)
print(dotplot(gsea_cp_sfg, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggplot2::labs(title = "GSEA Dotplot — Superior Frontal Gyrus"))
dev.off()

cat("GSEA plots saved ✅\n")
cat("\n=== 04_enrichment.R complete! ===\n")

