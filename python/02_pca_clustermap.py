# ============================================================
# 02_pca_clustermap.py
# Purpose: PCA, clustermap and DEG summary visualization
#          Comparison target: R (ggplot2 PCA, ComplexHeatmap)
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# --- Plot style ---
sns.set_theme(style="whitegrid", font_scale=1.2)

# --- Paths ---
FIGURE_DIR = "../figures/python_output/"
DATA_DIR   = "../data/processed/"

# ============================================================
# LOAD DATA
# ============================================================

top50 = pd.read_csv(DATA_DIR + "top50_matrix.csv", index_col=0)
meta  = pd.read_csv(DATA_DIR + "metadata.csv", index_col=0)
ec    = pd.read_csv(DATA_DIR + "deg_Entorhinal_Cortex.csv")
sfg   = pd.read_csv(DATA_DIR + "deg_Superior_Frontal_Gyrus.csv")

# Clean metadata — same pipeline as R and 01_eda_seaborn.py
meta['disease'] = meta['Disease State:ch1'].str.strip()
meta['region']  = meta['Organ Region:ch1'].str.strip()
meta['disease'] = meta['disease'].str.replace("Alzheimer's Disease", "AD")
meta_clean = meta.dropna(subset=['disease', 'region'])
meta_clean = meta_clean[meta_clean['disease'].isin(['AD', 'normal'])]

print("Data loaded ✅")
print(f"Clean samples: {len(meta_clean)}")

# ============================================================
# PLOT 1 — PCA (Seaborn + Matplotlib)
# ============================================================

# top50 shape: 50 genes × 105 samples
# For PCA we need samples as rows, genes as columns → transpose
# StandardScaler scales each gene to mean=0, std=1
# This is equivalent to scale. = TRUE in R's prcomp()

X = top50.T  # shape: 105 samples × 50 genes

# Scale the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Run PCA — n_components=2 gives us PC1 and PC2
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(X_scaled)

# Variance explained by each PC
var_explained = pca.explained_variance_ratio_ * 100

# Build a dataframe for plotting
# Merge PCA coordinates with metadata
pca_df = pd.DataFrame({
    'PC1':     pca_coords[:, 0],
    'PC2':     pca_coords[:, 1],
    'sample':  X.index
})

# Match sample order to metadata
pca_df = pca_df.set_index('sample')
pca_df = pca_df.join(meta_clean[['disease', 'region']])
pca_df = pca_df.dropna()

# Color by disease, marker by region
disease_colors = {'AD': '#E74C3C', 'normal': '#3498DB'}
region_markers = {
    'Entorhinal Cortex':     'o',
    'Hippocampus':           '^',
    'Medial Temporal Gyrus': 's',
    'Posterior Cingulate':   'P',
    'Primary Visual Cortex': 'X',
    'Superior Frontal Gyrus':'*'
}

fig, ax = plt.subplots(figsize=(10, 8))

for region, marker in region_markers.items():
    for disease, color in disease_colors.items():
        mask = (pca_df['region'] == region) & (pca_df['disease'] == disease)
        if mask.sum() > 0:
            ax.scatter(pca_df.loc[mask, 'PC1'],
                       pca_df.loc[mask, 'PC2'],
                       c=color,
                       marker=marker,
                       s=80,
                       alpha=0.8,
                       label=f"{region} — {disease}" if disease == 'AD' else "")

# Build clean legend
disease_patches = [mpatches.Patch(color=c, label=d)
                   for d, c in disease_colors.items()]
region_handles  = [plt.Line2D([0], [0], marker=m, color='grey',
                               linestyle='None', markersize=8, label=r)
                   for r, m in region_markers.items()
                   if r in pca_df['region'].values]

legend1 = ax.legend(handles=disease_patches,
                    title="Disease", loc='upper left', fontsize=9)
ax.add_artist(legend1)
ax.legend(handles=region_handles,
          title="Brain Region", loc='lower left', fontsize=8)

ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)", fontsize=12)
ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)", fontsize=12)
ax.set_title("PCA — Top 50 EC DEGs\nPython / Seaborn + scikit-learn",
             fontsize=13, fontweight='bold')
ax.axhline(0, color='grey', linewidth=0.5, linestyle='--')
ax.axvline(0, color='grey', linewidth=0.5, linestyle='--')

plt.tight_layout()
plt.savefig(FIGURE_DIR + "06_pca.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 06_pca.png ✅")

# ============================================================
# PLOT 2 — CLUSTERMAP (Seaborn)
# ============================================================

# Seaborn clustermap automatically clusters both rows and columns
# rows = genes, columns = samples
# Comparison to R's ComplexHeatmap

# Build color annotation for columns (samples)
# Map disease to color for the column color bar
disease_lut    = {'AD': '#E74C3C', 'normal': '#3498DB'}
col_colors     = meta_clean.loc[top50.columns, 'disease'].map(disease_lut)

# Build region color bar
region_palette = sns.color_palette("Set2", n_colors=6)
regions        = meta_clean['region'].unique()
region_lut     = dict(zip(regions, region_palette))
row_colors     = meta_clean.loc[top50.columns, 'region'].map(region_lut)

# Stack both annotation bars
col_colors_df  = pd.DataFrame({
    'Disease': col_colors,
    'Region':  row_colors
})

# Draw clustermap
# z_score=0 standardizes each row (gene) — same as scaling in R
cg = sns.clustermap(top50,
                    z_score=0,           # scale genes (rows)
                    cmap='coolwarm',     # blue=low, red=high
                    center=0,
                    col_colors=col_colors_df,
                    figsize=(16, 12),
                    xticklabels=False,   # too many samples to label
                    yticklabels=True,
                    dendrogram_ratio=0.1,
                    cbar_pos=(0.02, 0.8, 0.03, 0.15))

cg.ax_heatmap.set_title("Clustermap — Top 50 EC DEGs\nPython / Seaborn",
                          fontsize=13, fontweight='bold', pad=120)
cg.ax_heatmap.tick_params(axis='y', labelsize=8)

# Add legend for disease colors
disease_handles = [mpatches.Patch(color=c, label=d)
                   for d, c in disease_lut.items()]
cg.ax_col_dendrogram.legend(handles=disease_handles,
                             title="Disease",
                             loc='upper right',
                             fontsize=9)

plt.savefig(FIGURE_DIR + "07_clustermap.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 07_clustermap.png ✅")

# ============================================================
# PLOT 3 — REGION-WISE DEG SUMMARY (Seaborn barplot)
# ============================================================

# Count up/down DEGs per region and plot as grouped barplot
# Quick visual summary of the entire DE analysis

def count_degs(df, region_name):
    up   = ((df['adj.P.Val'] < 0.05) & (df['logFC'] >  1)).sum()
    down = ((df['adj.P.Val'] < 0.05) & (df['logFC'] < -1)).sum()
    return {'Region': region_name, 'Direction': 'Upregulated',   'Count': up}, \
           {'Region': region_name, 'Direction': 'Downregulated', 'Count': down}

rows = []
for record in count_degs(ec,  'Entorhinal Cortex'):
    rows.append(record)
for record in count_degs(sfg, 'Superior Frontal Gyrus'):
    rows.append(record)

summary_df = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(9, 6))

sns.barplot(data=summary_df,
            x='Region',
            y='Count',
            hue='Direction',
            palette={'Upregulated':   '#E74C3C',
                     'Downregulated': '#3498DB'},
            ax=ax)

# Add count labels on top of bars
for container in ax.containers:
    ax.bar_label(container, fontsize=10, padding=3)

ax.set_title("DEG Count per Brain Region — AD vs normal\nPython / Seaborn",
             fontsize=13, fontweight='bold')
ax.set_xlabel("Brain Region", fontsize=11)
ax.set_ylabel("Number of DEGs", fontsize=11)
ax.legend(title="Direction", fontsize=10)

plt.tight_layout()
plt.savefig(FIGURE_DIR + "08_deg_summary_barplot.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 08_deg_summary_barplot.png ✅")

print("\n=== 02_pca_clustermap.py complete! ===")
print("Check figures/python_output/ for plots 06-08")