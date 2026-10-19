import os
import numpy as np
import pandas as pd
import gseapy as gp
import scanpy as sc


np.random.seed(123)

# Dictionary to store variance ratios for each cell line
variance_data = {}


for cell_line in ['4T1','EO771']:

    # Get list of STING-dependent genes from interaction term in cGAMP vs H2O: STOP vs REF
    STING_df = pd.read_csv(''.join([cell_line,'_interaction/fcHurdle_STOP.csv']))
    STING_genes = STING_df.primerid[STING_df.fdr<0.5]
    
    # Read cGAMP adata
    adata = sc.read_h5ad(''.join(['adata_',cell_line,'_cGAMP.h5ad']))
    
    # Do PCs
    sc.pp.pca(adata)
    
    # We take the first 10 PCs
    variance_data[cell_line] = adata.uns['pca']['variance_ratio'][:10]
    
    #Create a DataFrame for PCA loadings (genes x PCs)
    n_PC = 10 # Number of PCs to save loadings for
    loadings_df = pd.DataFrame(
        adata.varm['PCs'], 
        index=adata.var_names, 
        columns=[f'PC{i+1}' for i in range(adata.varm['PCs'].shape[1])]
    ).iloc[:,0:n_PC]
    
    # Pick only STING-dependent genes
    loadings_df = loadings_df[loadings_df.index.isin(STING_genes)]
    
    PC_gsea_score = pd.DataFrame()
    
    for PC in ['PC1','PC2','PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']:
    
        # Run prerank
        ss = gp.prerank(
            rnk=loadings_df[PC].sort_values(ascending=False), 
            gene_sets='mouse_hallmarks.gmt', 
            outdir=None, 
            min_size=5,
            max_size=10000,
            permutation_num=1e4,
            seed=123, 
            sample_norm_method='custom', 
            no_plot=True
        )
        
        # Access the results (res2d typically contains the enrichment metrics)
        res = ss.res2d.copy()
        
        # Add the variant name to the results
        res['Name'] = PC
        
        # Select and append specific columns to the master dataframe
        PC_gsea_score = pd.concat([PC_gsea_score, res[['Name', 'Term', 'ES', 'NES',
                                                   'NOM p-val', 'FDR q-val',
                                                   'FWER p-val']]],
                                  ignore_index=True)
    PC_gsea_score.to_csv(''.join([cell_line,
                                  "_PC_cGAMP_gsea.csv"]), index=False)
    
# --- FINAL STEP: Save the Variance CSV ---
variance_df = pd.DataFrame(variance_data, index=[f'PC{i+1}' for i in range(10)])
variance_df.to_csv('PC_variance_explained.csv')