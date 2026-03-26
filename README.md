# 🧠 NeuroExprViz
### Cross-Tool Alzheimer's Transcriptomics Visualization

[![Live Dashboard](https://img.shields.io/badge/Live%20Dashboard-GitHub%20Pages-blue?style=flat-square&logo=github)](https://mayankgandhi13.github.io/NeuroExprViz/dashboard/index.html)
[![R](https://img.shields.io/badge/R-4.4.2-276DC3?style=flat-square&logo=r)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.11-3776AB?style=flat-square&logo=python)](https://www.python.org/)
[![JavaScript](https://img.shields.io/badge/JavaScript-ES6-F7DF1E?style=flat-square&logo=javascript)](https://developer.mozilla.org/en-US/docs/Web/JavaScript)
[![Dataset](https://img.shields.io/badge/Dataset-GSE5281-red?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)

> Differential expression analysis of Alzheimer's disease transcriptomics (GSE5281) — comparing visualization efficiency across **R**, **Python**, and **JavaScript** for bioinformatics workflows.

---

## 🔗 [→ Open Live Dashboard](https://mayankgandhi13.github.io/NeuroExprViz/dashboard/index.html)

---

## Project Overview

This project answers one question: **which visualization stack is most efficient for bioinformaticians?**

The same DE analysis results are reproduced across three ecosystems:

| Tool | Stack | Output |
|------|-------|--------|
| **R** | limma · ggplot2 · ComplexHeatmap · EnhancedVolcano | Publication-quality static figures |
| **Python** | Seaborn · Matplotlib · scikit-learn | EDA notebooks and comparison figures |
| **JavaScript** | Plotly.js | Interactive web dashboard |

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        DATA ACQUISITION                         │
│                                                                 │
│   NCBI GEO (GSE5281)  ──►  GEOquery (R)  ──►  ExpressionSet   │
│   161 samples · 6 brain regions · Affymetrix HG-U133 Plus 2.0  │
└─────────────────────────┬───────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                     PREPROCESSING & QC                          │
│                                                                 │
│   Raw Intensities  ──►  RMA Normalization  ──►  Log2 Matrix    │
│   54,675 probes    ──►  Probe Mapping      ──►  21,358 genes   │
│                         (hgu133plus2.db)                        │
│                    ──►  QC Plots (boxplot · density · PCA)     │
└─────────────────────────┬───────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│               DIFFERENTIAL EXPRESSION (limma)                   │
│                                                                 │
│   Region-stratified analysis · 5 brain regions                 │
│   Design: ~ 0 + disease  ──►  Contrast: AD - normal            │
│   eBayes shrinkage  ──►  BH correction  ──►  adj.P < 0.05      │
│                                                                 │
│   Entorhinal Cortex   →  3,358 DEGs (1,434↑ · 1,924↓)         │
│   Superior Frontal    →    841 DEGs   (691↑ ·   150↓)         │
│   Primary Visual      →      0 DEGs  (biological control)      │
└──────────────┬──────────────────────────┬───────────────────────┘
               │                          │
               ▼                          ▼
┌──────────────────────┐    ┌─────────────────────────────────────┐
│  ENRICHMENT ANALYSIS │    │           CSV EXPORTS               │
│                      │    │                                     │
│  GO ORA              │    │  deg_Entorhinal_Cortex.csv          │
│  clusterProfiler     │    │  deg_Superior_Frontal_Gyrus.csv     │
│                      │    │  gsea_EC.csv                        │
│  GSEA (Hallmark)     │    │  gsea_SFG.csv                       │
│  fgsea               │    │  top50_matrix.csv                   │
│                      │    │  metadata.csv                       │
└──────────────────────┘    └──────────┬──────────────────────────┘
                                       │
               ┌───────────────────────┼───────────────────────┐
               │                       │                       │
               ▼                       ▼                       ▼
┌──────────────────────┐ ┌─────────────────────┐ ┌────────────────────┐
│    R FIGURES         │ │  PYTHON FIGURES      │ │  JS DASHBOARD      │
│    (RStudio)         │ │  (VS Code)           │ │  (VS Code)         │
│                      │ │                      │ │                    │
│  EnhancedVolcano     │ │  Matplotlib volcano  │ │  Plotly.js volcano │
│  ComplexHeatmap      │ │  Seaborn clustermap  │ │  GSEA bubbles      │
│  GSEA ridge/dot      │ │  Seaborn PCA         │ │  DEG bar chart     │
│  PCA biplot          │ │  Correlation heatmap │ │  Gene search       │
│  ORA dotplot         │ │  logFC distribution  │ │  Region toggle     │
│                      │ │                      │ │                    │
│  14 PNGs             │ │  8 PNGs              │ │  GitHub Pages ✓    │
└──────────────────────┘ └─────────────────────┘ └────────────────────┘
```

---

## Dataset

| Property | Value |
|----------|-------|
| **Accession** | [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) |
| **Disease** | Alzheimer's Disease vs. Normal |
| **Samples** | 161 total (74 control · 87 AD) |
| **Platform** | Affymetrix HG-U133 Plus 2.0 (Microarray) |
| **Brain Regions** | Entorhinal Cortex · Hippocampus · Medial Temporal Gyrus · Posterior Cingulate · Primary Visual Cortex · Superior Frontal Gyrus |
| **Clean Samples** | 105 (after metadata filtering) |

---

## Key Biological Findings

**Entorhinal Cortex (most disrupted):**
- 3,358 DEGs — highest signal, first region affected in AD
- GO enrichment: synaptic vesicle cycle, vesicle-mediated transport
- GSEA: mitochondrial dysfunction (↓), neuroinflammation (↑)

**Superior Frontal Gyrus:**
- 841 DEGs — moderate disruption
- GO enrichment: TGF-beta signaling, glial cell differentiation
- GSEA: gliogenesis, astrocyte activation (↑)

**Primary Visual Cortex:**
- 0 DEGs — serves as biological negative control ✅

---

## Repository Structure

```
NeuroExprViz/
│
├── R/                          ← Open in RStudio
│   ├── 01_download_geo.R       # Download GSE5281 from NCBI GEO
│   ├── 02_preprocess_qc.R      # RMA normalization + QC plots
│   ├── 03_limma_de.R           # Region-stratified DE analysis
│   ├── 04_enrichment.R         # GO ORA + GSEA + enrichment plots
│   └── 05_figures.R            # Volcano · Heatmap · PCA
│
├── python/                     ← Open in VS Code
│   ├── 01_eda_seaborn.py       # Volcano · distribution · violin · heatmap
│   └── 02_pca_clustermap.py    # PCA · clustermap · DEG summary
│
├── dashboard/                  ← Deployed to GitHub Pages
│   ├── index.html              # Main dashboard page
│   ├── style.css               # Mixed theme styling
│   ├── volcano.js              # Interactive volcano + gene search
│   ├── summary.js              # DEG bar chart + GSEA bubble
│   └── comparison.js           # R vs Python vs JS tab comparison
│
├── data/
│   ├── raw/                    # GEO download cache (gitignored)
│   └── processed/              # CSV outputs from R pipeline
│
└── figures/
    ├── r_output/               # 14 PNGs from R
    └── python_output/          # 8 PNGs from Python
```

---

## Tool Comparison

| Feature | R | Python | JavaScript |
|---------|---|--------|------------|
| Setup complexity | Low | Low | Medium |
| Bioinformatics packages | **Excellent** | Good | None |
| Statistical analysis | **Excellent** | Good | None |
| Interactivity | None | Limited | **Full** |
| Web deployment | Shiny only | Dash/Streamlit | **Native** |
| ML integration | Good | **Excellent** | None |
| Publication figures | **Excellent** | Good | Limited |
| Lines of code (volcano) | ~8 | ~35 | ~80 |

**Conclusion:** R wins for analysis and publication figures. Python bridges bioinformatics and ML. JavaScript is essential for sharing interactive results with non-bioinformaticians.

---

## How to Run

**Prerequisites:** R 4.4+, Python 3.11+, modern browser

**R Pipeline:**
```bash
# Open R/ folder in RStudio
# Run scripts in order:
source("R/01_download_geo.R")
source("R/02_preprocess_qc.R")
source("R/03_limma_de.R")
source("R/04_enrichment.R")
source("R/05_figures.R")
```

**Python EDA:**
```bash
cd python/
conda activate binf6400
python 01_eda_seaborn.py
python 02_pca_clustermap.py
```

**Dashboard (local):**
```bash
# From project root
python -m http.server 8000
# Open: http://localhost:8000/dashboard/index.html
```

---

## Dependencies

**R packages:** GEOquery · Biobase · affy · limma · biomaRt · hgu133plus2.db · ggplot2 · ComplexHeatmap · EnhancedVolcano · clusterProfiler · fgsea · enrichplot · org.Hs.eg.db · msigdbr

**Python packages:** pandas · numpy · matplotlib · seaborn · scikit-learn · scipy

**JavaScript:** Plotly.js 2.27.0 (CDN)

---

## References

- Liang WS et al. (2008). Altered neuronal gene expression in brain regions differentially affected by Alzheimer's disease. *Physiol Genomics*. GSE5281.
- Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*.
- Wu T et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*.

---

<div align="center">
  <b>Mayank Gandhi</b> · Bioinformatics Portfolio<br>
  <a href="https://mayankgandhi13.github.io/NeuroExprViz/dashboard/index.html">Live Dashboard</a> ·
  <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281">Dataset (GSE5281)</a>
</div>
