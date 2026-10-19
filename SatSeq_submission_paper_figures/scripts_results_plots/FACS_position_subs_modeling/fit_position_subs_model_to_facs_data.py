import numpy as np
import pandas as pd
import sys

# Join path to functions file
sys.path.append('../common_use_functions')

from position_subs_effect_inference import *

#%% '''FACS STING, ATF3, IFN data'''
filename = '../../datasets/FACS/FACS_ATF3_IFN_NFkB_STING_Proliferation_4T1_EO771_all.xlsx'

with pd.ExcelFile(filename) as xls:
    for functional_readout in ['NFkB','ATF3','Dual_IFN','STING','Proliferation']:
        
        sheet_up = ''.join(['4T1_',functional_readout,'_FACS_All_vs_Up'])
        sheet_down = ''.join(['4T1_',functional_readout,'_FACS_All_vs_Down'])
        
        if functional_readout=='Proliferation':
            sheet_up = '4T1_Initial_vs_MitogenHigh_p20'
            sheet_down = '4T1_Initial_vs_MitogenLow_p20'        
        
        # Read sheets for up and down
        df_down = pd.read_excel(xls,sheet_name = sheet_down)
        df_up = pd.read_excel(xls,sheet_name = sheet_up)
        
        # Extract Replicate number
        patterns = {'Proliferation': (r"4T1_p20_(\d+)_", r"4T1_2FBS_p20_(\d+)_"),
                    'ATF3':          (r"ATF3_R(\d+)_Up", r"ATF3_R(\d+)_Dn"),
                    'Dual_IFN':      (r"R(\d+)_Up",      r"R(\d+)_Down"),
                    'STING':         (r"Up_(\d+)",      r"Down_(\d+)"),
                    'NFkB':          (r"NFkB_UP(\d+)",      r"NFkB_DN(\d+)")}

        pattern_up, pattern_down = patterns.get(functional_readout, (None, None))
            
        df_up["rep"] = df_up["experiment"].str.extract(pattern_up)[0]
        df_up["rep"] = pd.to_numeric(df_up["rep"]).fillna(1).astype(int)
        
        df_down["rep"] = df_down["experiment"].str.extract(pattern_down)[0]
        df_down["rep"] = pd.to_numeric(df_down["rep"]).fillna(1).astype(int)
        
        # --- Merge on variant and replicate ---
        merged = (
            df_up
            .merge(
                df_down,
                on=["variant", "rep", "is_control"],
                suffixes=("_Up", "_Down")
            )
        )
        
        # --- Build output dataframe ---
        fitness_data = pd.DataFrame({
            "variant": merged["variant"],
            "is_control": merged["is_control"],
            "fitness_Up": merged["observed_fitness_Up"],
            "fitness_Down": merged["observed_fitness_Down"],
            "experiment": (
                "Up_" + merged["rep"].astype(str) +
                "_Down_" + merged["rep"].astype(str)
            )
        })
        
        fitness_data['observed_fitness'] = fitness_data['fitness_Up'] - fitness_data['fitness_Down']
        
        # Use only df_up for proliferation
        if functional_readout=='Proliferation':
            fitness_data['observed_fitness'] = fitness_data['fitness_Up']            
        
        # Run fitting
        weights,aa_additive,aa_interaction,intercept_exp = position_subs_effect_inference(fitness_data,
                                                                       seed=123,
                                                                       advi_steps=80000)
        
        # Save inferred parameters
        weights.to_csv(''.join(['inferred_weights_FACS_4T1_',functional_readout,'.csv']))
        aa_additive.to_csv(''.join(['inferred_aamat_additive_FACS_4T1_',functional_readout,'.csv']))
        aa_interaction.to_csv(''.join(['inferred_aamat_interaction_FACS_4T1_',functional_readout,'.csv']))
        intercept_exp.to_csv(''.join(['inferred_intercept_experiment_FACS_4T1_',functional_readout,'.csv']))
    
    
    
    
    
