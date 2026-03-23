# ============================================================
# 03_limma_de.R
# Purpose: Clean metadata, filter samples, run region-
#          stratified differential expression with limma
# ============================================================

library(Biobase)
library(limma)
library(magrittr)
library(dplyr)
library(stringr)

# ============================================================
# JOB 1 — LOAD & CLEAN METADATA
# ============================================================

norm_expr <- readRDS("../data/processed/norm_expr_matrix.rds")
pdata     <- readRDS("../data/processed/metadata.rds")

# Extract the two columns we need
pdata$disease <- pdata$`Disease State:ch1`
pdata$region  <- pdata$`Organ Region:ch1`

# Fix region naming inconsistencies
# hippocampus → Hippocampus (capitalize)
# Posterior Singulate → Posterior Cingulate (typo fix)
pdata$region <- gsub("hippocampus", "Hippocampus", pdata$region)
pdata$region <- gsub("Posterior Singulate", "Posterior Cingulate", pdata$region)

# Fix disease label — remove apostrophe which causes R naming issues
pdata$disease <- gsub("Alzheimer's Disease", "AD", pdata$disease)

# Keep only samples with complete metadata
# nchar(trimws()) > 0 is more robust than != ""
# because it handles non-ASCII whitespace too
complete_idx <- !is.na(pdata$disease)           &
  !is.na(pdata$region)             &
  nchar(trimws(pdata$disease)) > 0 &
  nchar(trimws(pdata$region))  > 0

pdata_clean <- pdata[complete_idx, ]

# Use iconv + str_trim to strip non-ASCII whitespace characters
# that are commonly baked into GEO metadata files
pdata_clean$disease <- str_trim(iconv(pdata_clean$disease,
                                      from = "UTF-8",
                                      to   = "ASCII//TRANSLIT"))
pdata_clean$region  <- str_trim(iconv(pdata_clean$region,
                                      from = "UTF-8",
                                      to   = "ASCII//TRANSLIT"))

# Subset expression matrix to match the 105 clean samples
# rownames(pdata_clean) are the GEO sample IDs e.g. GSM119615
norm_expr_clean <- norm_expr[, rownames(pdata_clean)]

cat("=== After cleaning ===\n")
cat("Samples:", ncol(norm_expr_clean), "\n")
cat("Genes:", nrow(norm_expr_clean), "\n")
cat("\nSamples per region:\n")
print(table(pdata_clean$region))
cat("\nDisease labels:\n")
print(table(pdata_clean$disease))

# ============================================================
# JOB 2 — REGION-STRATIFIED LIMMA DE ANALYSIS
# ============================================================

regions     <- unique(pdata_clean$region)
deg_results <- list()

for (reg in regions) {
  
  cat("\n=== Processing:", reg, "===\n")
  
  # --- Step 1: Subset samples for this region ---
  region_idx <- pdata_clean$region == reg
  pdata_reg  <- pdata_clean[region_idx, ]
  expr_reg   <- norm_expr_clean[, rownames(pdata_reg)]
  
  cat("Samples in this region:", ncol(expr_reg), "\n")
  cat("AD:", sum(pdata_reg$disease == "AD"),
      "| normal:", sum(pdata_reg$disease == "normal"), "\n")
  
  # Skip region if it doesn't have both AD and normal samples
  # limma cannot estimate a contrast if one group is missing entirely
  if (sum(pdata_reg$disease == "AD") == 0 |
      sum(pdata_reg$disease == "normal") == 0) {
    cat("Skipping", reg, "— missing one group\n")
    next  # next skips to the next iteration of the loop
  }
  
  # --- Step 2: Build design matrix ---
  # levels = c("normal", "AD") ensures positive logFC = higher in AD
  # ~ 0 + disease means no intercept, one coefficient per group
  disease <- factor(pdata_reg$disease, levels = c("normal", "AD"))
  design  <- model.matrix(~ 0 + disease)
  colnames(design) <- levels(disease)
  
  # --- Step 3: Fit linear model ---
  # lmFit fits one linear model per gene across all samples
  fit <- lmFit(expr_reg, design)
  
  # --- Step 4: Define contrast AD vs normal ---
  # positive logFC = higher expression in AD vs normal
  contrasts <- makeContrasts(AD - normal, levels = design)
  fit2      <- contrasts.fit(fit, contrasts)
  
  # --- Step 5: Apply eBayes ---
  # borrows variance information across all genes
  # trend = TRUE accounts for mean-variance relationship in microarray
  fit2 <- eBayes(fit2, trend = TRUE)
  
  # --- Step 6: Extract results ---
  # number = Inf returns ALL genes not just top 10
  # BH = Benjamini-Hochberg FDR correction
  results      <- topTable(fit2,
                           coef          = 1,
                           number        = Inf,
                           adjust.method = "BH",
                           sort.by       = "P")
  results$gene <- rownames(results)
  
  # Count significant DEGs
  sig <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 1)
  cat("Significant DEGs (adj.P<0.05, |logFC|>1):", sig, "\n")
  
  # Store in list and save as CSV
  deg_results[[reg]] <- results
  filename <- paste0("../data/processed/deg_",
                     gsub(" ", "_", reg), ".csv")
  write.csv(results, filename, row.names = FALSE)
  cat("Saved:", filename, "✅\n")
}

cat("\n=== All regions processed! ===\n")

# ============================================================
# JOB 3 — SUMMARY TABLE
# ============================================================

cat("\n=== DEG Summary per Region ===\n")
for (reg in regions) {
  
  # Skip regions that were not processed
  if (is.null(deg_results[[reg]])) {
    cat(reg, "→ Skipped (missing group)\n")
    next
  }
  
  res  <- deg_results[[reg]]
  up   <- sum(res$adj.P.Val < 0.05 & res$logFC >  1)
  down <- sum(res$adj.P.Val < 0.05 & res$logFC < -1)
  cat(reg, "→ Up:", up, "| Down:", down, "| Total:", up + down, "\n")
}

