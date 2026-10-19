import numpy as np
import pandas as pd
import sys

# Join path to functions file
sys.path.append('../common_use_functions')

from position_subs_effect_inference import *

# Define datasets, suffixes, and whether they require standard deviation calculation
datasets = [
    {"file": "Mammary_Fat_Pad_Primary_vs_Met_boot_4T1.csv", "suffix": "Primary_vs_Met_4T1", "return_sd": True},
    {"file": "Mammary_Fat_Pad_Initial_vs_Met_boot_4T1.csv", "suffix": "Initial_vs_Met_4T1", "return_sd": True},
    {"file": "Mammary_Fat_Pad_Initial_vs_Primary_boot_4T1.csv", "suffix": "Initial_vs_Primary_4T1", "return_sd": False},
    {"file": "Mammary_Fat_Pad_Initial_vs_Primary_boot_EO771.csv", "suffix": "Initial_vs_Primary_EO771", "return_sd": False},
    {"file": "Tail_Vein_vs_Met_boot_4T1.csv", "suffix": "Tail_Vein_vs_Met_4T1", "return_sd": True},
    {"file": "Initial_vs_Tail_Vein_boot_EO771.csv", "suffix": "Initial_vs_Tail_Vein_EO771", "return_sd": False}
]

# Loop over the files
for ds in datasets:
    print(f"Processing: {ds['suffix']}...")
    
    # Load dataset
    fitness_data = pd.read_csv(f"../../datasets/bulk_fitness/{ds['file']}")
    
    if ds["return_sd"]:
        # Unpack 5 returned values when return_sd=True
        weights, weights_sd, aa_additive, aa_interaction, intercept_exp = position_subs_effect_inference(
            fitness_data, seed=123, advi_steps=30000, return_sd=True
        )
        # Save the extra standard deviation file
        weights_sd.to_csv(f"inferred_weights_sd_{ds['suffix']}.csv")
    else:
        # Unpack 4 returned values when return_sd is default/False
        weights, aa_additive, aa_interaction, intercept_exp = position_subs_effect_inference(
            fitness_data, seed=123, advi_steps=30000
        )
        
    # Save standard inferred parameters common to all runs
    weights.to_csv(f"inferred_weights_{ds['suffix']}.csv")
    aa_additive.to_csv(f"inferred_aamat_additive_{ds['suffix']}.csv")
    aa_interaction.to_csv(f"inferred_aamat_interaction_{ds['suffix']}.csv")
    intercept_exp.to_csv(f"inferred_intercept_experiment_{ds['suffix']}.csv")
