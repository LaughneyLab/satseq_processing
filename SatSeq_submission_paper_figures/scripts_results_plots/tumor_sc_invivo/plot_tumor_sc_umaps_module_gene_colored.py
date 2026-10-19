import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap

# ==============================================================================
# 1. LOAD DATA (Assuming sequential cell ordering matches perfectly)
# ==============================================================================
# Read single-cell object
adata = sc.read_h5ad("adata_4T1_tumor.h5ad")

# Get STING-dependent genes
sting_df = pd.concat([
    pd.read_csv('4T1_interaction/fcHurdle_STOP.csv'),
    pd.read_csv('EO771_interaction/fcHurdle_STOP.csv')
])

# 2. Extract unique primerids where FDR < 0.5
sting_genes = sting_df.loc[sting_df['fdr'] < 0.5, 'primerid'].unique()

# Use only STING genes
# Create a shallow copy and subset AnnData to save memory overhead
adata_sub = adata[:, sting_genes].copy()
adata_sub = adata_sub[adata_sub.obs['variant_substitution'].to_numpy() != 'None']

adata = adata_sub

# Read module scores CSV (saved with .reset_index(drop=True))
module_scores = pd.read_csv('hotspot_module_scores_tumor_4T1.csv', index_col=0)

# ==============================================================================
# 2. RUN EMBEDDING AND ALIGN POSITIONALLY
# ==============================================================================
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Slice adata purely by length to ensure a strict 1-to-1 sequential match with the CSV 
adata_subset = adata[:len(module_scores)].copy()
umap_coords = adata_subset.obsm['X_umap']

# ==============================================================================
# 3. EXTRACT MODULE SCORES SECURELY BY POSITION
# ==============================================================================
cols = [c for c in module_scores.columns if c != 'variant_substitution']

col_prolif = cols[0]  
col_ifn    = cols[1]  
col_nfkb   = cols[2]  

# Bypass label-based index matching entirely
m_prolif = module_scores[col_prolif].values
m_ifn    = module_scores[col_ifn].values
m_nfkb   = module_scores[col_nfkb].values

# ==============================================================================
# 3. DEFINE LAYOUT & COLORMAPS
# ==============================================================================
# 5 Columns represent: Col 0 = Hotspot Score, Cols 1-4 = Specific Marker Genes
fig = plt.figure(figsize=(35, 12))

# Left margin is shifted inward (left=0.08) to safely hold the vertical colorbar
gs = fig.add_gridspec(
    3, 5, 
    hspace=0.35, 
    wspace=0.1, 
    top=0.90, 
    bottom=0.05,
    left=0.08,
    right=0.98
)

# Custom grayscale colormap for all gene expression panels
gray_cmap = LinearSegmentedColormap.from_list('custom_gray',
                                              ['#E5E4E2', '#708090'])

# Define matrix structure of axes for easy loop access
axes = np.empty((3, 5), dtype=object)
for r in range(3):
    for c in range(5):
        axes[r, c] = fig.add_subplot(gs[r, c])

# ==============================================================================
# 4. PLOT ELEVATED COLORBAR (Hotspot module score)
# ==============================================================================
norm_fixed = mcolors.Normalize(vmin=-10.0, vmax=10.0)
# Position the colorbar vertically on the far left side of the canvas
# Coordinates: [left, bottom, width, height]
ax_cb = fig.add_axes([0.02, 0.1, 0.012, 0.3]) 

cb_rdbu = ColorbarBase(ax_cb, cmap='RdBu_r', norm=norm_fixed, orientation='vertical')
ax_cb.set_title("Hotspot Module score", fontsize=25, weight='bold',
                pad=20, rotation=90)
ax_cb.yaxis.set_label_position('left') # Force text to the left of the colorbar ticks
ax_cb.tick_params(labelsize=25)

# ==============================================================================
# 5. PLOT COLUMN 0: HOTSPOT SCORES
# ==============================================================================
# Proliferation Score
axes[0, 0].scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_prolif, cmap='RdBu_r', vmin=-10, vmax=10, s=2, edgecolors='none')
axes[0, 0].set_title("Proliferation", fontsize=25, weight='bold', pad=15)

# IFN Score
axes[1, 0].scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_ifn, cmap='RdBu_r', vmin=-10, vmax=10, s=2, edgecolors='none')
axes[1, 0].set_title("IFN", fontsize=25, weight='bold', pad=15)

# NFkB Score
axes[2, 0].scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_nfkb, cmap='RdBu_r', vmin=-10, vmax=10, s=2, edgecolors='none')
axes[2, 0].set_title("NFkB", fontsize=25, weight='bold', pad=15)

# Clear boundaries for Column 0
for r in range(3):
    axes[r, 0].set_axis_off()

# ==============================================================================
# 6. PLOT COLUMNS 1-4: GENE EXPRESSION GRID (Grayscale & Italics)
# ==============================================================================
# Define the gene matrix exactly mapping to Rows (0-2) and Columns (1-4)
gene_grid = [
    ["Mki67",  "Cdc20",  "Cdc6",   "Cks1b"], # Row 0: Proliferation markers
    ["Cxcl10", "Isg15",  "Ifit3",  "Irf1" ], # Row 1: IFN markers
    ["Nfkb1",  "Atf3",   "Nfkbia", "Rel"  ]  # Row 2: NFkB markers
]

# Populate gene UMAPs using Scanpy engine
for r in range(3):
    for c in range(1, 5):
        target_gene = gene_grid[r][c-1]
        
        # Plot directly into the designated grid cell axis
        sc.pl.umap(
            adata_subset, 
            color=target_gene, 
            size=2, 
            frameon=False, 
            show=False, 
            ax=axes[r, c], 
            cmap=gray_cmap, 
            title=""
        )
        
        # Override title aesthetics: No bolding, styled in raw text Italics
        axes[r, c].set_title(target_gene, fontsize=25, style='italic', pad=15)

# ==============================================================================
# 7. EXPORT IMAGE FILE
# ==============================================================================
plt.savefig('hotspot_module_gene_umaps_tumor_4T1.png', dpi=200, bbox_inches='tight')
plt.show()