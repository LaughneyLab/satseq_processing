import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

# 1. Read single-cell data
adata = sc.read_h5ad("../mt269/satseq/matt_cbsulaughney/10x/lung_mets/adata_save_2026_07_05.h5ad")

# 2. Extract base names from the 'sample' column safely
adata.obs["variant"] = adata.obs["sample"].astype(str).str.split("_").str.get(0)

# 3. Rename specific cell types in the 'celltype' column
celltype_map = {
    "AM": "Alveolar macrophages",
    "Endothelial cells": "Endothelial",
    "IM": "Interstitial macrophages",
    "ILC": "Innate lymphoid cells",
    "iMON": "Inflammatory monocytes",
    "SCMF": "Secondary crest myofibroblasts"
}
adata.obs["celltype"] = adata.obs["celltype"].astype(str).replace(celltype_map).astype("category")

# 4. Define the custom hex color palette for variants
variant_colors = {
    "P13": "#FBB040",
    "K337": "#E43B3F",
    "REF": "#4E4E4E",
    "STOP": "#B4B4B4"
}

# 5. Enforce categorical ordering and assign variant colors
adata.obs["variant"] = adata.obs["variant"].astype("category")
categories = adata.obs["variant"].cat.categories
color_list = [variant_colors.get(cat, "#D3D3D3") for cat in categories]
adata.uns["variant_colors"] = color_list

# 6. Initialize layout
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 7))

# 7. Plot Left UMAP (Variant)
sc.pl.umap(
    adata, 
    color="variant", 
    frameon=False, 
    show=False, 
    save=False,
    legend_loc=None,  
    ax=axes[0]
)
axes[0].set_title("") 
axes[0].set_aspect('equal', adjustable='box')

# Manually build custom bottom legend handles for Variant
variant_cats = adata.obs["variant"].cat.categories
variant_palette = adata.uns["variant_colors"]
var_handles = [
    Line2D([], [], marker='o', color='w', markerfacecolor=color, markersize=8, label=cat)
    for cat, color in zip(variant_cats, variant_palette)
]
axes[0].legend(
    handles=var_handles,
    loc="upper center", 
    bbox_to_anchor=(0.5, -0.05), 
    ncol=4,  
    frameon=False,
    fontsize=9
)

# 8. Plot Right UMAP (Celltype) without built-in on-data text (to manually control it)
sc.pl.umap(
    adata, 
    color="celltype", 
    frameon=False, 
    show=False, 
    save=False,
    legend_loc=None,  # Disabled built-in text to prevent standard overlapping
    ax=axes[1]
)
axes[1].set_title("") 
axes[1].set_aspect('equal', adjustable='box')

# 9. Manually calculate cluster centers and apply overlap repulsion
texts = []
umap_coords = adata.obsm["X_umap"]

for celltype in adata.obs["celltype"].cat.categories:
    # Identify indices belonging to the current cell type
    mask = adata.obs["celltype"] == celltype
    if not np.any(mask):
        continue
    
    # Calculate median X and Y coordinates for the cluster center
    x_pos, y_pos = np.median(umap_coords[mask], axis=0)
    
    # Create the text object on the plot canvas
    t = axes[1].text(
        x_pos, y_pos, celltype, 
        fontsize=8, 
        weight="normal",
        bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1) # Semi-transparent background
    )
    texts.append(t)

# Automatically push text boxes apart and draw subtle connecting line pointers
adjust_text(
    texts, 
    ax=axes[1], 
    arrowprops=dict(arrowstyle="-", color="black", lw=0.5, alpha=0.6),
    expand=(1.2, 1.4)  # Padding multipliers to push names slightly further apart if crowded
)

# 10. Fine-tune layout parameters
plt.subplots_adjust(wspace=0.05, bottom=0.15, left=0.05, right=0.95)

# 11. Export high-res file
fig.savefig("umap_lung_mets_variant_colored.png", dpi=300, bbox_inches="tight")
plt.close(fig)
