import numpy as np
import pandas as pd
import arviz as az
import pymc as pm


filename = '../../datasets/FACS/FACS_ATF3_IFN_NFkB_STING_Proliferation_4T1_EO771_all.xlsx'

with pd.ExcelFile(filename) as xls:
    for functional_readout in ['STING']:
        
        sheet_up = ''.join(['4T1_',functional_readout,'_FACS_All_vs_Up'])
        sheet_down = ''.join(['4T1_',functional_readout,'_FACS_All_vs_Down'])
        
        # Read sheets for up and down
        df_down = pd.read_excel(xls,sheet_name = sheet_down)
        df_up = pd.read_excel(xls,sheet_name = sheet_up)
        
        # Extract Replicate number
        patterns = {'STING':         (r"Up_(\d+)",      r"Down_(\d+)")}

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

        
# Parse variant col to get WT and mutant AA and mutation position
fitness_data[['wt_residue', 'position', 'mutant_residue']] = fitness_data[
    'variant'].str.extract(r'^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$')

fitness_data['position'] = pd.to_numeric(fitness_data['position'])

fitness_data[['wt_residue', 'mutant_residue']] = fitness_data[
    ['wt_residue', 'mutant_residue']].astype(str)

# Drop, e.g., '_379X','dCTT' residues
fitness_data = fitness_data[fitness_data[['wt_residue',
                                          'position',
                                          'mutant_residue']].isna().any(axis=1)==False]

#AA name list
aa_list=np.array(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
                 'S','T','V','W','Y'])

fitness_data = fitness_data[fitness_data[['wt_residue',
                                          'mutant_residue']].isin(aa_list).all(axis=1)==True]

# Drop unmutated residues, if any
fitness_data = fitness_data[fitness_data.wt_residue!=fitness_data.mutant_residue]

#%% '''Featurize and index mouse and mutations'''

mutation_position =  fitness_data.position.to_numpy().astype(int) - 1 

# List which AA substitutions are in data (based on full AA matrix)
# Indices of index and mutant AAs in AA list (flattened)
wt_aa = np.where(np.array(list(''.join(fitness_data['wt_residue'])))[:, None] == 
                    aa_list[None, :])[1]

mutant_aa = np.where(np.array(list(''.join(fitness_data['mutant_residue'])))[:, None] == 
                    aa_list[None, :])[1]

aa_subs_index = (20*wt_aa + mutant_aa)
aa_subs_index[wt_aa==mutant_aa]=0 #assign 0 to positions with unchanged AA

# Renumber AA change index matrix with unique AA change indices
unique_aa_subs,_, aa_sub_index_unique = np.unique(aa_subs_index, 
                               return_index=True,
                               return_inverse=True);

# Index mouse_experiment combos
experiment_list,_, experiment_index = np.unique(fitness_data.experiment, 
                               return_index=True,
                               return_inverse=True)

# Read fitness data
mutant_fitness = fitness_data['observed_fitness']
#%% '''Set up hierarchical model'''
with pm.Model() as mutant_fitness_model: 
    
    # position-indepedent common amino acid distance matrix multiplier flattened
    aa_subs_effect = pm.Normal("aa_effect",
                            mu = pm.Normal("mu_aa", mu=0, sigma=1),
                            sigma = pm.Exponential("sigma_aa",lam=1),
                            shape=(len(unique_aa_subs),2))    
                      
    # positional weights, one for each mouse/experiment
    weights = pm.Normal("weights",
                      mu = pm.Normal("mu_w", mu=0, sigma=1),
                      sigma = pm.Exponential("sigma_w",lam=1),
                      shape = (378,2))
    
    # Experiment_specific intercept
    intercept_experiment = pm.Normal("intercept_experiment",
                            mu = pm.Normal("mu_exp", mu=0, sigma=1),
                            sigma = pm.Exponential("sigma_exp",lam=1),
                            shape = (len(experiment_list)))
    
    # Predictor
    predicted_fitness = intercept_experiment[experiment_index] + weights[
        mutation_position,0] + aa_subs_effect[aa_sub_index_unique,0] + weights[
            mutation_position,1] * aa_subs_effect[aa_sub_index_unique,1] 
    
    # Fit
    fitness_obs = pm.Normal("fitness_obs",
                       mu=predicted_fitness,
                       sigma=pm.Exponential("sigma_obs",lam=1),
                       observed=mutant_fitness) 
    
# Generate the plate diagram using ArviZ
graph = pm.model_to_graphviz(mutant_fitness_model) 
graph.render("mutant_fitness_model", format="png", cleanup=True)
