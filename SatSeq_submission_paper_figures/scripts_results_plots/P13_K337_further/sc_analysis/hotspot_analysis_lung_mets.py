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
from matplotlib.lines import Line2D


'''Process sc data'''

# Get STING-dependent genes
sting_df = pd.concat([
    pd.read_csv('fcHurdle_4T1_interaction_STOP.csv'),
    pd.read_csv('fcHurdle_EO771_interaction_STOP.csv')
])

sting_genes = sting_df.loc[sting_df['fdr'] < 0.5, 'primerid'].unique()

df_cgamp_stop = pd.read_csv('fcHurdle_4T1_cGAMP_STOP.csv')
filtered_genes = df_cgamp_stop.loc[df_cgamp_stop['coef'] < 0, 'primerid'].unique()
sting_genes = list(set(sting_genes) & set(filtered_genes))

# Read sc data
adata = sc.read_h5ad("../mt269/satseq/matt_cbsulaughney/10x/lung_mets/adata_save_2026_07_05.h5ad")

# Subset for tumor sc
adata = adata[adata.obs["celltype"] == "4T1", :].copy()

# Use only STING genes
# Create a shallow copy and subset AnnData to save memory overhead
adata_sub = adata[:, sting_genes].copy()

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
lcz = hs.compute_local_correlations(hs_genes, jobs=4)

modules = hs.create_modules(
    min_gene_threshold=32, core_only=True, fdr_threshold=0.05)

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
hs.results.join(hs.modules).to_csv('hotspot_gene_modules_lung_mets.csv')

# 4. Reorder rows and columns of the matrix to match the heatmap
ordered_z_matrix = hs.local_correlation_z.loc[gene_order, gene_order]

# 5. Export the structured matrix to a CSV file
ordered_z_matrix.to_csv('hotspot_gene_distances_lung_mets.csv')

# Aggregate module score per cell
module_scores = hs.calculate_module_scores()
module_scores['sample'] = adata.obs['sample']
module_scores.reset_index(drop=True).to_csv('hotspot_module_scores_4T1_lung_mets.csv')

# Plot module scores on umap
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# 1. Align adata with module_scores
common_cells = module_scores.index
adata_subset = adata[common_cells].copy()
umap_coords = adata_subset.obsm['X_umap']

# 2. Dynamically extract module columns
module_cols = [c for c in module_scores.columns if c != 'sample']
num_modules = len(module_cols)
num_panels = num_modules + 1  # Total panels = All modules + 1 Variant panel

# 3. Extract and map variants
adata_subset.obs["variant"] = adata_subset.obs["sample"].astype(str).str.split("_").str.get(0)
variant_colors = {
    "P13": "#FBB040",
    "K337": "#E43B3F",
    "REF": "#4E4E4E",
    "STOP": "#B4B4B4"
}
var_colors_mapped = [variant_colors.get(v, "#000000") for v in adata_subset.obs["variant"]]

# Point size configuration
pt_size = 0.8

# 4. Set up layout dynamically based on panel count
# Dynamically scale figure width (4.5 inches per panel matches original proportions)
fig_width = 4.5 * num_panels
fig = plt.figure(figsize=(fig_width, 7.5))
gs = fig.add_gridspec(2, num_panels, height_ratios=[5, 0.8], hspace=0.35, wspace=0.18, bottom=0.15)

# --- 5. LOOP AND PLOT MODULES ---
for i, col_name in enumerate(module_cols):
    ax = fig.add_subplot(gs[0, i])
    
    # Extract scores and fill missing values
    scores = np.nan_to_num(module_scores.loc[adata_subset.obs_names, col_name].values, nan=-5.0)
    
    # Plot module UMAP
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=scores, cmap='RdBu_r', vmin=-5, vmax=5, s=pt_size, edgecolors='none')
    ax.set_title(f"Module {i + 1}", fontsize=13, weight='bold', pad=8)
    ax.set_axis_off()

# --- 6. PLOT VARIANT (Rightmost Panel) ---
variant_idx = num_panels - 1
ax_var = fig.add_subplot(gs[0, variant_idx])
ax_var.scatter(umap_coords[:, 0], umap_coords[:, 1], c=var_colors_mapped, s=pt_size, edgecolors='none')
ax_var.set_title("Variant", fontsize=13, weight='bold', pad=8)
ax_var.set_axis_off()

# --- 7. SCALE BARS & LEGENDS (Row 1) --- 
norm_fixed = mcolors.Normalize(vmin=-5.0, vmax=5.0)

# Shared horizontal module colorbar placed under the approximate middle module panel
middle_idx = num_modules // 2
ax_rdbu_cb = fig.add_subplot(gs[1, middle_idx])
cb_rdbu = ColorbarBase(ax_rdbu_cb, cmap='RdBu_r', norm=norm_fixed, orientation='horizontal')
cb_rdbu.set_label("Module Score", fontsize=11, labelpad=8)
cb_rdbu.ax.tick_params(labelsize=10)

# Categorical Legend underneath Variant panel (Rightmost column)
ax_leg = fig.add_subplot(gs[1, variant_idx])
ax_leg.set_axis_off()

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='P13', markerfacecolor=variant_colors['P13'], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='K337', markerfacecolor=variant_colors['K337'], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='STOP', markerfacecolor=variant_colors['STOP'], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='REF', markerfacecolor=variant_colors['REF'], markersize=8)
]

ax_leg.legend(handles=legend_elements, loc='center', ncol=4, frameon=False, fontsize=11, handletextpad=0.2, columnspacing=1.0)

plt.savefig('hotspot_module_score_umaps_lung_mets.png', dpi=200, bbox_inches='tight')
plt.show()


# Do EnrichR
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
            organism="mouse",
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
output_excel = "hotspot_modules_gsea_lung_mets.xlsx"

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    wide_pvalue.to_excel(writer, sheet_name="Adjusted P-value")
    wide_overlap.to_excel(writer, sheet_name="Overlap")
    wide_genes.to_excel(writer, sheet_name="Genes")

    







