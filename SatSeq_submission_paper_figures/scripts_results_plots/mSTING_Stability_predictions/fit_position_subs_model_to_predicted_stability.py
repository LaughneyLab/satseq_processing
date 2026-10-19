import numpy as np
import pandas as pd
import sys

# Join path to functions file
sys.path.append('../common_use_functions')

from position_subs_effect_inference import *

# 1. Define the WT sequence
wt_seq = "".join(["MPYSNLHPAIPRPRGHRSKYVALIFLVASLMILWVAKDPPNHTLKYLALHLASHELGLLL",
                 "KNLCCLAEELCHVQSRYQGSYWKAVRACLGCPIHCMAMILLSSYFYFLQNTADIYLSWMF",
                 "GLLVLYKSLSMLLGLQSLTPAEVSAVCEEKKLNVAHGLAWSYYIGYLRLILPGLQARIRM",
                 "FNQLHNNMLSGAGSRRLYILFPLDCGVPDNLSVVDPNIRFRDMLPQQNIDRAGIKNRVYS",
                 "NSVYEILENGQPAGVCILEYATPLQTLFAMSQDAKAGFSREDRLEQAKLFCRTLEEILED",
                 "VPESRNNCRLIVYQEPTDGNSFSLSQEVLRHIRQEEKEEVTMNAPMTSVAPPPSVLSQEP",
                 "RLLISGMDQPLPLRTDLI"])

# Extract ddG data
def process_and_plot_heatmap(filename,
                             pos_col,
                             mutant_col,
                             ddg_col,
                             scale=1.0):
    # 2. Read and clean data
    df = pd.read_csv(filename)
    
    df_clean = pd.DataFrame()
    df_clean['position'] = df[pos_col].astype(int)
    df_clean['mutant_residue'] = df[mutant_col].astype(str)
    df_clean['observed_fitness'] = df[ddg_col].astype(float) * scale
    
    # Filter positions
    df_clean = df_clean[(df_clean['position'] >= 1) & (df_clean['position'] <= 378)]
    
    # Extract WT character and create 'variant' column (WT-position-mutant)
    # Pandas 0-index requires position - 1
    df_clean['wt_residue'] = df_clean['position'].apply(lambda p: wt_seq[p - 1])
    df_clean['variant'] = (df_clean['wt_residue'] + 
                           df_clean['position'].astype(str) + 
                           df_clean['mutant_residue'])
    df_clean['experiment'] = "ddG_model"
    return df_clean


model_configs = {
    "ThermoMPNN": {
        "filename": "ThermoMPNN/output_mSTING_AF_dimer_ThermoMPNN_Single_threshold_100_Distance_5_Penalize.csv",
        "model_name": "ThermoMPNN",
        "pos_col": "Position",
        "mutant_col": "Mutant AA",
        "ddg_col": "ddG (kcal/mol)",
        "scale": -1.0
    },
    "ProteinMPNN": {
        "filename": "ProteinMPNN-ddG/output_mSTING_AF_dimer_ProteinMPNN-ddG_model_v_48_020_rep_5_seed_123.csv",
        "model_name": "ProteinMPNN",
        "pos_col": "Position",
        "mutant_col": "Post",
        "ddg_col": "ProteinMPNN-ddG Score",
        "scale": 1.0  # Default scale explicitly added
    },
    "DeepDDG": {
        "filename": "DeepDDG/output_cols_split.csv",
        "model_name": "DeepDDG",
        "pos_col": "ResID",
        "mutant_col": "Mut",
        "ddg_col": "ddG",
        "scale": 1.0  # Default scale explicitly added
    }
}

# Loop over models
for model in ["ThermoMPNN","ProteinMPNN","DeepDDG"]:
    
    # Select the specific model configuration you want to run
    config = model_configs[model]

    # Call the function by passing the dictionary values explicitly
    df_ddg = process_and_plot_heatmap(
        filename=config["filename"],
        pos_col=config["pos_col"],
        mutant_col=config["mutant_col"],
        ddg_col=config["ddg_col"],
        scale=config["scale"])
    
    # Run positional weights model
    weights, aa_additive, aa_interaction, intercept_exp = position_subs_effect_inference(
        df_ddg, seed=123, advi_steps=30000)
    
    # Save standard inferred parameters common to all runs
    weights.to_csv(f"inferred_weights_{model}.csv")
    aa_additive.to_csv(f"inferred_aamat_additive_{model}.csv")
    aa_interaction.to_csv(f"inferred_aamat_interaction_{model}.csv")
    
    




