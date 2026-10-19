# Follows https://hotspot.readthedocs.io/en/latest/CD4_Tutorial.html
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.utils import sparsefuncs
from scipy.sparse import issparse
import re
import hotspot
import gseapy as gp
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase

'''Process sc data'''

# Get STING-dependent genes
sting_df = pd.concat([
    pd.read_csv('4T1_interaction/fcHurdle_STOP.csv'),
    pd.read_csv('EO771_interaction/fcHurdle_STOP.csv')
])

# 2. Extract unique primerids where FDR < 0.5
sting_genes = sting_df.loc[sting_df['fdr'] < 0.5, 'primerid'].unique()

# Read sc data
adata = sc.read_h5ad("adata_4T1_tumor.h5ad")

# Use only STING genes
# Create a shallow copy and subset AnnData to save memory overhead
adata_sub = adata[:, sting_genes].copy()
adata_sub = adata_sub[adata_sub.obs['variant_substitution'].to_numpy() != 'None']

# Filter zero variance genes
mean, gene_vars = sparsefuncs.mean_variance_axis(adata_sub.X, axis=0)
adata_sub = adata_sub[:, gene_vars > 0]

# --- OPTIMIZED TRANSFORMATION ---
# 1. Access the sparse matrix directly
raw_counts = adata_sub.X.copy()

# 2. Inverse log transform non-zero entries while keeping matrix sparse
if issparse(raw_counts):
    raw_counts.data = np.expm1(raw_counts.data)
else:
    raw_counts = np.expm1(raw_counts)

# 3. Scale by library size
scaling_factor = adata_sub.obs['library_size'].values / 1e4
raw_counts = raw_counts.multiply(scaling_factor[:, np.newaxis])

# 4. Finalize: Round, convert, and store
raw_counts.data = np.rint(raw_counts.data)
adata_sub.layers["counts_csc"] = raw_counts.tocsc()

'''Run Hotspot and extract per cell module scores'''

hs = hotspot.Hotspot(
    adata_sub, 
    layer_key="counts_csc", 
    model='danb', 
    latent_obsm_key="X_pca",
    umi_counts_obs_key="library_size"
)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30,
)

# compute autocorrelations for each gene, in the pca-space
hs_results = hs.compute_autocorrelations(jobs=1)

# Select the genes with significant lineage autocorrelation
hs_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values('Z', ascending=False).head(500).index

# Compute pair-wise local correlations between these genes
# Sometimes need to run this in console separately if jobs>1
lcz = hs.compute_local_correlations(hs_genes, jobs=1)

modules = hs.create_modules(
    min_gene_threshold=64, core_only=True, fdr_threshold=0.05)

# Plot gene modules
hs.plot_local_correlations(vmin=-100, vmax=100)

# Extract gene order in heatmap
# 1. Get the original input gene labels used for clustering
genes_input = hs.local_correlation_z.columns

# 2. Extract the clustered index order from Hotspot's linkage tree
# (no_plot=True prevents an extra matplotlib plot from rendering)
dendro = dendrogram(hs.linkage, no_plot=True)
clustered_indices = dendro['leaves']

# 3. Reorder the gene names to get the exact heatmap order
gene_order = genes_input[clustered_indices].tolist()


# Show top genes in a module
for module in np.arange(1,np.max(modules)+1):
    print(module)
    results = hs.results.join(hs.modules)
    results = results.loc[results.Module == module]
    print(results.sort_values('Z', ascending=False).index[0:30])

# Save module genes
hs.results.join(hs.modules).to_csv('hotspot_gene_modules.csv')

# 4. Reorder rows and columns of the matrix to match the heatmap
ordered_z_matrix = hs.local_correlation_z.loc[gene_order, gene_order]

# 5. Export the structured matrix to a CSV file
ordered_z_matrix.to_csv('hotspot_gene_distances.csv')

# Aggregate module score per cell
module_scores = hs.calculate_module_scores()
module_scores['variant_substitution'] = adata.obs['variant_substitution']
module_scores.reset_index(drop=True).to_csv('hotspot_module_scores_tumor_4T1.csv')

# Plot module scores on umap
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, frameon=False)

''' Plot UMAPS '''
# 1. Align adata with module_scores
common_cells = module_scores.index
adata_subset = adata[common_cells].copy()
umap_coords = adata_subset.obsm['X_umap']

# 2. Extract columns strictly by position index
cols = [c for c in module_scores.columns if c != 'variant_substitution']

col_prolif = cols[0]  
col_ifn    = cols[1]  
col_nfkb   = cols[2]  

m_prolif = module_scores.loc[adata_subset.obs_names, col_prolif].values
m_ifn    = module_scores.loc[adata_subset.obs_names, col_ifn].values
m_nfkb   = module_scores.loc[adata_subset.obs_names, col_nfkb].values

m_prolif = np.nan_to_num(m_prolif, nan=-5.0)
m_ifn    = np.nan_to_num(m_ifn, nan=-5.0)
m_nfkb   = np.nan_to_num(m_nfkb, nan=-5.0)

# Helper function to scale raw scores from [-5, 5] cleanly for blending
def scale_to_range(x):
    return np.clip((x - (-5.0)) / (5.0 - (-5.0)), 0, 1)

w_prolif = scale_to_range(m_prolif)
w_ifn    = scale_to_range(m_ifn)
w_nfkb   = scale_to_range(m_nfkb)

# 3. Define the anchor colors for the 4th Blended Panel
hex_base     = "#e0e0e0"  
hex_charcoal = "#6082B6"  
hex_red      = "#e3312a"  
hex_green    = "#7ac943"  

rgb_base     = np.array(mcolors.to_rgb(hex_base))
rgb_charcoal = np.array(mcolors.to_rgb(hex_charcoal))
rgb_red      = np.array(mcolors.to_rgb(hex_red))
rgb_green    = np.array(mcolors.to_rgb(hex_green))

# 4. Mix colors linearly for the 4th Blended Panel
blended_rgb = rgb_base + (
    w_prolif[:, None] * (rgb_charcoal - rgb_base) +
    w_ifn[:, None] * (rgb_red - rgb_base) +
    w_nfkb[:, None] * (rgb_green - rgb_base)
)
blended_rgb = np.clip(blended_rgb, 0.0, 1.0)

# 5. Set up layout with Gridspec (Increased height and margins for larger text)
fig = plt.figure(figsize=(18, 6.5))
gs = fig.add_gridspec(2, 4, height_ratios=[5, 0.7], hspace=0.3, wspace=0.15)

# Main UMAP Axis panels (Row 0)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

# Panel 1: Proliferation (Font size increased)
ax1.scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_prolif, cmap='RdBu_r', vmin=-5, vmax=5, s=2, edgecolors='none')
ax1.set_title(f"Proliferation\n({col_prolif})", fontsize=13, weight='bold', pad=8)
ax1.set_axis_off()

# Panel 2: IFN (Font size increased)
ax2.scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_ifn, cmap='RdBu_r', vmin=-5, vmax=5, s=2, edgecolors='none')
ax2.set_title(f"IFN\n({col_ifn})", fontsize=13, weight='bold', pad=8)
ax2.set_axis_off()

# Panel 3: NFkB (Font size increased)
ax3.scatter(umap_coords[:, 0], umap_coords[:, 1], c=m_nfkb, cmap='RdBu_r', vmin=-5, vmax=5, s=2, edgecolors='none')
ax3.set_title(f"NFkB\n({col_nfkb})", fontsize=13, weight='bold', pad=8)
ax3.set_axis_off()

# Panel 4: Blended Mix (Font size increased)
ax4.scatter(umap_coords[:, 0], umap_coords[:, 1], c=blended_rgb, s=2, edgecolors='none')
ax4.set_title("Blended Triplet\n(All Mixed)", fontsize=13, weight='bold', pad=8)
ax4.set_axis_off()

# --- 6. ADDING THE SCALE BARS (Row 1 with larger fonts) ---
norm_fixed = mcolors.Normalize(vmin=-5.0, vmax=5.0)

# A. Shared horizontal colorbar under Panel 2 representing the 'RdBu_r' scale
ax_rdbu_cb = fig.add_subplot(gs[1, 1])  
cb_rdbu = ColorbarBase(ax_rdbu_cb, cmap='RdBu_r', norm=norm_fixed, orientation='horizontal')
cb_rdbu.set_label("Module Score (Panels 1-3)", fontsize=11, labelpad=8)
cb_rdbu.ax.tick_params(labelsize=10)

# B. 3 stacked horizontal scale bars underneath Panel 4 using a sub-gridspec
gs_blended_cbs = gs[1, 3].subgridspec(3, 1, hspace=0.6)
ax_cb1 = fig.add_subplot(gs_blended_cbs[0, 0])
ax_cb2 = fig.add_subplot(gs_blended_cbs[1, 0])
ax_cb3 = fig.add_subplot(gs_blended_cbs[2, 0])

cmap_prolif = mcolors.LinearSegmentedColormap.from_list("prolif_map", [hex_base, hex_charcoal])
cmap_ifn    = mcolors.LinearSegmentedColormap.from_list("ifn_map", [hex_base, hex_red])
cmap_nfkb   = mcolors.LinearSegmentedColormap.from_list("nfkb_map", [hex_base, hex_green])

cb1 = ColorbarBase(ax_cb1, cmap=cmap_prolif, norm=norm_fixed, orientation='horizontal')
cb1.set_ticks([]) 
ax_cb1.set_title("Proliferation", loc='left', fontsize=10, pad=2)

cb2 = ColorbarBase(ax_cb2, cmap=cmap_ifn, norm=norm_fixed, orientation='horizontal')
cb2.set_ticks([])
ax_cb2.set_title("IFN", loc='left', fontsize=10, pad=2)

cb3 = ColorbarBase(ax_cb3, cmap=cmap_nfkb, norm=norm_fixed, orientation='horizontal')
cb3.ax.tick_params(labelsize=9) # Clean [-5, 0, 5] tick sizing
ax_cb3.set_title("NFkB", loc='left', fontsize=10, pad=2)

plt.savefig('hotspot_module_score_umaps.png', dpi=200, bbox_inches='tight')

plt.show()

'''Do EnrichR'''

# 1. Master list to collect results from each comparison
master_results = []
results = hs.results.join(hs.modules) 

# 2. Iterate through each CSV file/module
for module in np.arange(1, np.max(modules) + 1):
    results_module = results.loc[results.Module == module]
    
    # Extract gene list (using the index since it contains the gene names)
    gene_list = results_module["Z"].sort_values(ascending=False).index.tolist()
    
    # 3. Skip if the list is empty to prevent Enrichr errors
    if not gene_list:
        continue

    # 4. Run Enrichr
    try:
        # Note: You can also use a custom .gmt file like in your prerank code
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets="mouse_hallmarks.gmt",
            organism="Mouse",
            outdir=None
        )
        
        # Extract the DataFrame results
        res = enr.results.copy()
        if res.empty:
            continue
            
        # Add tracking column for the module
        res["module"] = module 
        
        # Append to our master list
        master_results.append(res)
        
    except Exception as e:
        print(f"Error running Enrichr for module {module}: {e}")

# 5. Combine everything into a single dataframe
final_df = pd.concat(master_results, ignore_index=True)
final_df = final_df[['module', 'Term', 'Overlap', 'P-value', 'Adjusted P-value', 'Combined Score', 'Genes']]

# Create wide-form dataframes for your metrics of interest
wide_pvalue = final_df.pivot(index="Term", columns="module", values="Adjusted P-value").fillna(1)
wide_overlap = final_df.pivot(index="Term", columns="module", values="Overlap").fillna("0/0")
wide_genes = final_df.pivot(index="Term", columns="module", values="Genes").fillna("")

# Save all sheets into a single Excel file
output_excel = "hotspot_modules_gsea_tumor.xlsx"

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    wide_pvalue.to_excel(writer, sheet_name="Adjusted P-value")
    wide_overlap.to_excel(writer, sheet_name="Overlap")
    wide_genes.to_excel(writer, sheet_name="Genes")

    







