import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex, to_rgb
from adjustText import adjust_text

# 1. Set global random seeds for absolute consistency across runs
np.random.seed(42)
sc.settings.seed = 42

# 2. Read single-cell data
adata = sc.read_h5ad("../mt269/satseq/matt_cbsulaughney/10x/lung_mets/adata_save_2026_07_05.h5ad")

# 3. Extract base names from 'sample' column safely
adata.obs["variant"] = adata.obs["sample"].astype(str).str.split("_").str.get(0)

# 4. Create combined annotations column with T-cell subtype prioritization
celltype_series = adata.obs["celltype"].astype(str).str.strip()

if "T_cell_subtype" in adata.obs.columns:
    t_cell_mask = adata.obs["T_cell_subtype"].notna() & (adata.obs["T_cell_subtype"].astype(str) != "nan")
    celltype_series[t_cell_mask] = adata.obs["T_cell_subtype"].astype(str)[t_cell_mask].str.strip()

# 5. Apply case-insensitive cell type mapping
celltype_map = {
    "am": "Alveolar macrophages",
    "endothelial cells": "Endothelial",
    "endothelial": "Endothelial",
    "im": "Interstitial macrophages",
    "ilc": "Innate lymphoid cells",
    "imon": "Inflammatory monocytes",
    "scmf": "Secondary crest myofibroblasts",
    "4t1": "Tumor cells"
}

mapped_series = celltype_series.str.lower().map(celltype_map)
adata.obs["detailed_celltype"] = mapped_series.fillna(celltype_series).astype("category")

# 6. Define custom color palette for variant variables
variant_colors = {
    "P13": "#FBB040",
    "K337": "#E43B3F",
    "REF": "#4E4E4E",
    "STOP": "#B4B4B4"
}
adata.obs["variant"] = adata.obs["variant"].astype("category")
var_categories = adata.obs["variant"].cat.categories
adata.uns["variant_colors"] = [variant_colors.get(cat, "#D3D3D3") for cat in var_categories]

# 7. Generate a deterministic, vibrant color palette without yellow
ct_categories = adata.obs["detailed_celltype"].cat.categories
num_ct = len(ct_categories)

from scanpy.plotting import palettes
vibrant_pool = list(palettes.default_28) if hasattr(palettes, 'default_28') else list(palettes.default_20)

backup_vibrants = ["#1f77b4", "#9467bd", "#e377c2", "#17becf", "#2ca02c", "#bcbd22", "#d62728"]
for c in backup_vibrants:
    if c not in vibrant_pool:
        vibrant_pool.append(c)

final_ct_colors = []
for idx, cat in enumerate(ct_categories):
    if "4t1" in str(cat).lower():
        final_ct_colors.append("#4B0082")  # Indigo / Dark Purple
        continue
        
    color = vibrant_pool[idx % len(vibrant_pool)]  # Indexing dynamically fixes variance
    rgb = to_rgb(color)
    
    # Check individual RGB channels: high red, high green, low blue
    is_yellow = (rgb[0] > 0.70 and rgb[1] > 0.70 and rgb[2] < 0.40)
    if is_yellow:
        color = "#FF1493"  # Deep Pink
        
    final_ct_colors.append(to_hex(color))

adata.uns["detailed_celltype_colors"] = final_ct_colors


# =========================================================================
# PLOT 1: HIGH-RES PNG WITH CELL TYPE LABELS OVERLAID
# =========================================================================
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))

# Left Plot (Variant)
sc.pl.umap(adata, color="variant", frameon=False, show=False, save=False, legend_loc=None, ax=axes[0])
axes[0].set_title("")
axes[0].set_aspect('equal', adjustable='box')

var_handles = [
    Line2D([], [], marker='o', color='w', markerfacecolor=color, markersize=8, label=cat) 
    for cat, color in zip(var_categories, adata.uns["variant_colors"])
]
axes[0].legend(handles=var_handles, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, fontsize=9)

# Right Plot (Celltype)
sc.pl.umap(adata, color="detailed_celltype", frameon=False, show=False, save=False, legend_loc=None, ax=axes[1])
axes[1].set_title("")
axes[1].set_aspect('equal', adjustable='box')

# Calculate cluster centers and apply non-overlapping labels
texts = []
umap_coords = adata.obsm["X_umap"]
for celltype in ct_categories:
    mask = adata.obs["detailed_celltype"] == celltype
    if not np.any(mask): 
        continue
    x_pos, y_pos = np.median(umap_coords[mask], axis=0)
    t = axes[1].text(
        x_pos, y_pos, celltype, 
        fontsize=8, 
        weight="normal", 
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1)
    )
    texts.append(t)

adjust_text(texts, ax=axes[1], arrowprops=dict(arrowstyle="-", color="black", lw=0.5, alpha=0.6), expand=(1.3, 1.5))

plt.subplots_adjust(wspace=0.05, bottom=0.15, left=0.05, right=0.95)
fig.savefig("umap_lung_mets_variant_colored_T_subsets_included.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# =========================================================================
# PLOT 2: HIGH-RES CLEAN PNG WITHOUT ANY TEXT LABELS
# =========================================================================
fig_clean, axes_clean = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))

# Left Clean Plot
sc.pl.umap(adata, color="variant", frameon=False, show=False, save=False, legend_loc=None, ax=axes_clean[0])
axes_clean[0].set_title("")
axes_clean[0].set_aspect('equal', adjustable='box')
axes_clean[0].legend(handles=var_handles, loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, fontsize=9)

# Right Clean Plot
sc.pl.umap(adata, color="detailed_celltype", frameon=False, show=False, save=False, legend_loc=None, ax=axes_clean[1])
axes_clean[1].set_title("")
axes_clean[1].set_aspect('equal', adjustable='box')

plt.subplots_adjust(wspace=0.05, bottom=0.15, left=0.05, right=0.95)
fig_clean.savefig("umap_lung_mets_clean_nolabels.png", dpi=300, bbox_inches="tight")
plt.close(fig_clean)
