import os
# Get the path of the current module
current_module_path = os.path.dirname(__file__)
# Define the path to your data directory
data_directory = current_module_path

'''
Function to run hierarchical Bayesian model for position and aa subs effect inference
'''
def position_subs_effect_inference(fitness_data,
                                   seed=123,
                                   advi_steps=30000,
                                   return_sd = False):


    import numpy as np
    import pandas as pd
    import arviz as az
    import pymc as pm
    
    np.random.seed(seed)
    
    # Parse variant col to get WT and mutant AA and mutation position
    fitness_data[['wt_residue', 'position', 'mutant_residue']] = fitness_data[
        'variant'].str.extract(r'^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$')

    fitness_data['position'] = pd.to_numeric(fitness_data['position'])

    fitness_data[['wt_residue', 'mutant_residue']] = fitness_data[
        ['wt_residue', 'mutant_residue']].astype(str)

    # Drop, e.g., '_378X','dCTT' residues
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
        
       # Sampling with approximate posterior
        with mutant_fitness_model:
             posterior_draws=pm.fit(n=advi_steps,method="advi",
                                    random_seed=seed,progressbar=True)
             inferred_params = az.summary(posterior_draws.sample(10000,
                                                                 random_seed=seed))
        
    #%% '''Extract parameters and check fit'''
    # Find start index of inferred weight parameters
    weight_index = np.argwhere(inferred_params.index=='weights[0, 0]')[0,0]
        
    # extract and reshape positional weights
    inferred_weights=inferred_params.iloc[weight_index:weight_index+378*2,0
                                          ].to_numpy().reshape((378,2), order='C')
    
    # Make a df and add colnames
    inferred_weights = pd.DataFrame(inferred_weights,
                                    columns=['additive','interaction'])
    
    # SD of weights
    if return_sd:
        # extract and reshape positional weights
        inferred_weights_sd=inferred_params.iloc[weight_index:weight_index+378*2,1
                                              ].to_numpy().reshape((378,2), order='C')
        
        # Make a df and add colnames
        inferred_weights_sd = pd.DataFrame(inferred_weights_sd,
                                        columns=['additive','interaction'])

    # Find start index of inferred aa subs parameters
    aa_effect_index = np.argwhere(inferred_params.index=='aa_effect[0, 0]')[0,0]
        
    # extract and reshape aa_distance
    inferred_aa_effect=inferred_params.iloc[
        aa_effect_index:aa_effect_index+len(unique_aa_subs)*2,
        0].to_numpy().reshape((len(unique_aa_subs),2), order='C')
    
    # # Initialize with inferred mean of the pooling distribution for AA multiplier
    aa_sub_matrix_additive = np.zeros((20, 20)) + inferred_params.loc['mu_aa','mean']
    aa_sub_matrix_multiplicative = np.zeros((20, 20)) + inferred_params.loc['mu_aa','mean']

    for wt in range(20):
        for mut in range(20):
            
            if wt==mut:
                aa_sub_matrix_additive[wt, mut] = 0
                aa_sub_matrix_multiplicative[wt, mut] = 0
            
            if 20*wt+mut in unique_aa_subs:
                # Takes care of AA reindexing done at the beginning
                aa_sub_matrix_additive[wt, mut] = inferred_aa_effect[np.where(unique_aa_subs==20*wt+mut)[0][0],0]
                aa_sub_matrix_multiplicative[wt, mut] = inferred_aa_effect[np.where(unique_aa_subs==20*wt+mut)[0][0],1]

    inferred_aa_matrix_additive = pd.DataFrame(aa_sub_matrix_additive,
                                 index=aa_list,
                                 columns=aa_list)
                                     
    inferred_aa_matrix_interaction = pd.DataFrame(aa_sub_matrix_multiplicative,
                                 index=aa_list,
                                 columns=aa_list)

    # Find start index of inferred experiment parameters
    experiment_intercept_index = np.argwhere(inferred_params.index=='intercept_experiment[0]')[0,0]
        
    # extract inferred experiment parameters
    inferred_intercept_experiment=inferred_params.iloc[
        experiment_intercept_index:experiment_intercept_index+len(experiment_list),
        0].to_numpy()
    
    # Make into a df
    inferred_intercept_experiment = pd.DataFrame(inferred_intercept_experiment,
                                                 index=experiment_list)    
    
    if return_sd:
        return inferred_weights,inferred_weights_sd,inferred_aa_matrix_additive, inferred_aa_matrix_interaction,inferred_intercept_experiment
    else:
        return inferred_weights,inferred_aa_matrix_additive, inferred_aa_matrix_interaction,inferred_intercept_experiment
            
    
    
