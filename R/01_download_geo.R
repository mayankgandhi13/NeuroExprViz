# ============================================================
# 01_download_geo.R
# Purpose: Download GSE5281 from NCBI GEO and save locally
# ============================================================

library(GEOquery)
library(Biobase)

# --- Set paths ---
raw_dir <- "../data/raw"
dir.create(raw_dir, showWarnings = FALSE)

# --- Download GSE5281 ---
cat("Downloading GSE5281 from NCBI GEO...\n")

gse <- getGEO("GSE5281", 
              destdir = raw_dir,
              GSEMatrix = TRUE, 
              AnnotGPL = FALSE)

# GSE5281 comes as a list — extract the first element
eset <- gse[[1]]

# --- Quick sanity check ---
cat("=== Dataset Overview ===\n")
cat("Samples:", ncol(eset), "\n")
cat("Probes:", nrow(eset), "\n")
cat("Columns in metadata:", colnames(pData(eset)), "\n")

# --- Save locally so we never re-download ---
saveRDS(eset, file = "../data/raw/GSE5281_eset.rds")
cat("Saved to data/raw/GSE5281_eset.rds ✅\n")


