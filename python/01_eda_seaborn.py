# ============================================================
# 01_eda_seaborn.py
# Purpose: Reproduce and extend DEG visualization using Python
#          Comparison target: R (ggplot2, EnhancedVolcano)
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# --- Plot style ---
sns.set_theme(style="whitegrid", palette="muted", font_scale=1.2)

# --- Global color map — used across all plots ---
colors = {'Upregulated':   '#E74C3C',
          'Downregulated': '#3498DB',
          'NS':            '#95A5A6'}

# --- Paths ---
FIGURE_DIR = "../figures/python_output/"
DATA_DIR   = "../data/processed/"

# ============================================================
# LOAD DATA
# ============================================================

ec    = pd.read_csv(DATA_DIR + "deg_Entorhinal_Cortex.csv")
sfg   = pd.read_csv(DATA_DIR + "deg_Superior_Frontal_Gyrus.csv")
top50 = pd.read_csv(DATA_DIR + "top50_matrix.csv", index_col=0)
meta  = pd.read_csv(DATA_DIR + "metadata.csv", index_col=0)

print("Data loaded ✅")
print(f"EC genes: {len(ec)}")
print(f"SFG genes: {len(sfg)}")
print(f"Top50 matrix shape: {top50.shape}")

# ============================================================
# CATEGORIZE DEGs — shared helper used across plots
# ============================================================

def categorize(row):
    if row['adj.P.Val'] < 0.05 and row['logFC'] > 1:
        return 'Upregulated'
    elif row['adj.P.Val'] < 0.05 and row['logFC'] < -1:
        return 'Downregulated'
    else:
        return 'NS'

ec['category']  = ec.apply(categorize, axis=1)
sfg['category'] = sfg.apply(categorize, axis=1)

# ============================================================
# PLOT 1 — VOLCANO PLOT (Matplotlib)
# ============================================================

def plot_volcano(df, title, filename):

    fig, ax = plt.subplots(figsize=(10, 8))

    for cat, color in colors.items():
        mask = df['category'] == cat
        ax.scatter(df.loc[mask, 'logFC'],
                   -np.log10(df.loc[mask, 'adj.P.Val']),
                   c=color,
                   alpha=0.5,
                   s=15,
                   label=f"{cat} (n={mask.sum()})")

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(1,  color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(-1, color='grey', linestyle='--', linewidth=0.8)

    top_genes = df[df['category'] != 'NS'].nsmallest(15, 'adj.P.Val')
    for _, row in top_genes.iterrows():
        ax.annotate(row['gene'],
                    xy=(row['logFC'], -np.log10(row['adj.P.Val'])),
                    fontsize=7,
                    alpha=0.8)

    ax.set_xlabel("Log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-Log₁₀ (adj. P-value)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig(FIGURE_DIR + filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename} ✅")

plot_volcano(ec,
             "Volcano Plot — Entorhinal Cortex (AD vs normal)\nPython / Matplotlib",
             "01_volcano_EC.png")

plot_volcano(sfg,
             "Volcano Plot — Superior Frontal Gyrus (AD vs normal)\nPython / Matplotlib",
             "02_volcano_SFG.png")

# ============================================================
# PLOT 2 — logFC DISTRIBUTION (Seaborn)
# ============================================================

ec['region']  = 'Entorhinal Cortex'
sfg['region'] = 'Superior Frontal Gyrus'
combined = pd.concat([ec, sfg], ignore_index=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, (region, group) in zip(axes, combined.groupby('region')):
    sns.histplot(data=group,
                 x='logFC',
                 hue='category',
                 palette=colors,
                 kde=True,
                 bins=80,
                 ax=ax,
                 alpha=0.6)
    ax.set_title(region, fontsize=12, fontweight='bold')
    ax.set_xlabel("Log₂ Fold Change")
    ax.axvline(0,  color='black', linestyle='-',  linewidth=0.8)
    ax.axvline(1,  color='grey',  linestyle='--', linewidth=0.8)
    ax.axvline(-1, color='grey',  linestyle='--', linewidth=0.8)

plt.suptitle("logFC Distribution — AD vs normal", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(FIGURE_DIR + "03_logfc_distribution.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 03_logfc_distribution.png ✅")

# ============================================================
# PLOT 3 — VIOLIN PLOTS OF KEY AD GENES (Seaborn)
# ============================================================

norm = pd.read_csv(DATA_DIR + "norm_expr_matrix.csv", index_col=0)

# Clean metadata
meta['disease'] = meta['Disease State:ch1'].str.strip()
meta['region']  = meta['Organ Region:ch1'].str.strip()
meta['disease'] = meta['disease'].str.replace("Alzheimer's Disease", "AD")

meta_clean = meta.dropna(subset=['disease', 'region'])
meta_clean = meta_clean[meta_clean['disease'].isin(['AD', 'normal'])]

# Key AD biomarker genes
ad_genes = ['GFAP', 'NEAT1', 'APOE', 'CAMK1G', 'GABRD',
            'SST', 'NDUFA1', 'CD163', 'ABCA1', 'VIM']
ad_genes = [g for g in ad_genes if g in norm.index]
print(f"\nAD genes found in matrix: {ad_genes}")

# Subset and reshape
expr_subset = norm.loc[ad_genes, meta_clean.index].T
expr_subset = expr_subset.merge(meta_clean[['disease']],
                                left_index=True,
                                right_index=True)
expr_long = expr_subset.melt(id_vars='disease',
                              var_name='gene',
                              value_name='expression')

fig, ax = plt.subplots(figsize=(14, 6))

sns.violinplot(data=expr_long,
               x='gene',
               y='expression',
               hue='disease',
               palette={'AD': '#E74C3C', 'normal': '#3498DB'},
               split=True,
               inner='box',
               ax=ax)

ax.set_title("Expression of Key AD Genes — AD vs normal\nPython / Seaborn",
             fontsize=13, fontweight='bold')
ax.set_xlabel("Gene", fontsize=11)
ax.set_ylabel("Normalized Expression (log₂)", fontsize=11)
ax.tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(FIGURE_DIR + "04_violin_ad_genes.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 04_violin_ad_genes.png ✅")

# ============================================================
# PLOT 4 — CORRELATION HEATMAP (Seaborn)
# ============================================================

corr_matrix = top50.T.corr()

fig, ax = plt.subplots(figsize=(14, 12))

sns.heatmap(corr_matrix,
            cmap='coolwarm',
            center=0,
            square=True,
            linewidths=0.3,
            cbar_kws={'label': 'Pearson Correlation'},
            xticklabels=True,
            yticklabels=True,
            ax=ax)

ax.set_title("Gene-Gene Correlation — Top 50 EC DEGs\nPython / Seaborn",
             fontsize=13, fontweight='bold')
ax.tick_params(axis='x', rotation=90, labelsize=7)
ax.tick_params(axis='y', rotation=0,  labelsize=7)

plt.tight_layout()
plt.savefig(FIGURE_DIR + "05_correlation_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: 05_correlation_heatmap.png ✅")

print("\n=== 01_eda_seaborn.py complete! ===")
print("Check figures/python_output/ for 5 PNG files")