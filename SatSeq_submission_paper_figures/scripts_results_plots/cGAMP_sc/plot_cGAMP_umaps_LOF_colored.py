import scanpy as sc
import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

for cell_line in ["4T1","EO771"]:
    
    # Load sc data
    adata = sc.read_h5ad("".join(["adata_",cell_line,"_cGAMP.h5ad"]))
    
    # Filter out ambiguous variants
    adata = adata[adata.obs['variant_substitution'].to_numpy() != 'None']
    
    # Load LOF variants
    lof_df = pd.read_csv('LOF_variants.csv')
    lof_variants = lof_df['variant'].to_list()
    
    # 2. Extract base variant (e.g., P13X/P13A -> P13)
    # Extracts the letters and all subsequent digits (e.g., A141F -> A141, M1D -> M1)
    adata.obs['base_variant'] = adata.obs['variant_substitution'].str.extract(r'([A-Za-z]+\d+)')
    
    # 2. Handle special cases
    # Define conditions
    is_wt = adata.obs['variant_substitution'] == 'mSTING_WT'
    is_special = adata.obs['variant_substitution'].isin(['dCTT', 'STOP', 'RRDD'])
    
    # Apply logic: If WT -> REF; If special -> keep original name; Else -> keep regex result
    adata.obs['base_variant'] = np.select(
        [is_wt, is_special], 
        ['REF', adata.obs['variant_substitution']], 
        default=adata.obs['base_variant']
    )
    
    # 1. Create a column specifically for these four groups
    # We'll label everything else as 'Other' to keep the focus clear
    conditions = [
        adata.obs['base_variant'] == 'REF',
        #adata.obs['base_variant'] == 'dCTT',
        #adata.obs['base_variant'] == 'RRDD',
        adata.obs['base_variant'] == 'STOP'
        
    ]
    choices = ['REF',
               #'dCTT',
               #'RRDD',
               'STOP']
    
    adata.obs['controls'] = np.select(conditions, choices, default='Other')
    
    # 2. Define colors (using bright colors for the 4 groups and light gray for 'Other')
    custom_palette = {
        'REF': '#4d9221',    
        'STOP': '#c51b7d',   
        'Other': '#A9A9A9'   
    }
    
    # 1. Identify which variants belong to which LOF category
    ifna_list = lof_df[lof_df['group'].str.contains("IFNa", na=False)]['variant'].tolist()
    upr_list = lof_df[lof_df['group'].str.contains("UPR", na=False)]['variant'].tolist()
    nfkb_list = lof_df[lof_df['group'].str.contains("NFkB", na=False)]['variant'].tolist()
    
    # 2. Create indicator columns in adata (True/False as strings for coloring)
    adata.obs['IFN_LOF'] = adata.obs['base_variant'].isin(ifna_list).astype(str)
    adata.obs['UPR_LOF'] = adata.obs['base_variant'].isin(upr_list).astype(str)
    adata.obs['NFkB_LOF'] = adata.obs['base_variant'].isin(nfkb_list).astype(str)
    
    # 3. Define distinct palettes for each plot
    # Format: { 'True': Color, 'False': Background }
    ifn_pal = {'True': '#e3312a', 'False': '#ececec'}  # Red
    upr_pal = {'True': '#35b6e3', 'False': '#ececec'}  # Blue
    nfkb_pal = {'True': '#7ac943', 'False': '#ececec'} # Green
    
    # 3. Plot the 3 UMAPs
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    for ax in axes.ravel():
        ax.title.set_fontsize(25)       # Adjust this number to make UMAP titles larger or smaller

    fig.suptitle(f"{cell_line}", fontsize=25, x=0, y=0.98)
    
    # Top Row: LOF UMAPs
    sc.pl.umap(adata, color='controls', groups=['REF', 'STOP'], palette=custom_palette, 
               size=20, title='Controls', frameon=False, show=False, ax=axes[0, 0], na_in_legend=False)
    sc.pl.umap(adata, color='IFN_LOF', groups=['True'], palette=ifn_pal, 
               title='IFN LOF', size=20, frameon=False, show=False, ax=axes[0, 1], legend_loc=None)
    sc.pl.umap(adata, color='UPR_LOF', groups=['True'], palette=upr_pal, 
               title='UPR LOF', size=20, frameon=False, show=False, ax=axes[0, 2], legend_loc=None)
    sc.pl.umap(adata, color='NFkB_LOF', groups=['True'], palette=nfkb_pal, 
               title='NFkB LOF', size=20, frameon=False, show=False, ax=axes[0, 3], legend_loc=None)
    
    
    # Bottom Row: Gene Expression UMAPs
    # Create a custom grayscale colormap ranging from your custom light gray to black
    gray_cmap = LinearSegmentedColormap.from_list('custom_gray',
                                                  ['#D3D3D3', '#737373'])
    
    sc.pl.umap(adata, color='Tmem173', size=20, frameon=False, show=False,
               ax=axes[1, 0], cmap=gray_cmap, title="Sting1")
    sc.pl.umap(adata, color='Cxcl10', size=20, frameon=False, show=False,
               ax=axes[1, 1], cmap=gray_cmap)
    sc.pl.umap(adata, color='Atf3', size=20, frameon=False, show=False,
               ax=axes[1, 2], cmap=gray_cmap)
    sc.pl.umap(adata, color='Nfkb1', size=20, frameon=False, show=False,
               ax=axes[1, 3], cmap=gray_cmap)
    
    active_plots = [
    axes[0, 0], axes[0, 1], axes[0, 2], axes[0, 3], # Top row
    axes[1, 0], axes[1, 1], axes[1, 2], axes[1, 3]              # Bottom row
    ]   

    for ax in active_plots:
        # This explicitly overrides any title text Scanpy locked down
        current_text = ax.get_title()
        ax.set_title(current_text, fontsize=25, pad=15)
        
    leg = axes[0, 0].get_legend()
    if leg is not None:
        # Scale up individual legend entries ('REF', 'STOP')
        for text in leg.get_texts():
            text.set_fontsize(18)
    
    # Scale up the legend title header ('controls')
    leg.get_title().set_fontsize(20)
    
        
    plt.tight_layout()
    
    # Save in high resolution
    plt.savefig("".join(["umaps_",cell_line,"_cGAMP_LOF_colored.png"]),
                format='png', dpi=300, bbox_inches='tight')
    