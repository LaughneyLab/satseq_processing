def generate_single_mutant_dict(protein_seq):
    """
    Given a protein sequence (string), return a dictionary of all
    single amino-acid mutants named like 'G12D'. Includes WT.
    """
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")   # 20 standard amino acids
    seq = protein_seq.strip()
    L = len(seq)

    mutant_dict = {}

    # Add WT sequence
    mutant_dict["WT"] = seq

    # Loop over all positions
    for pos in range(L):
        wt_aa = seq[pos]
        for mut_aa in aa_list:
            if mut_aa == wt_aa:
                continue

            # Generate mutant sequence
            new_seq = list(seq)
            new_seq[pos] = mut_aa
            mutated_seq = "".join(new_seq)

            # Mutation name: WTaa + (1-based position) + MUTaa
            mut_name = f"{wt_aa}{pos+1}{mut_aa}"

            # Store in dictionary
            mutant_dict[mut_name] = mutated_seq

    return mutant_dict

# Get predictions in batches of 5
def generate_groups(start, end, group_size):
    numbers = list(range(start, end))
    for i in range(0, len(numbers), group_size):
        yield numbers[i:i + group_size]

import numpy as np
import pandas as pd
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()

# Mouse WT STING sequence from https://rest.uniprot.org/uniprotkb/Q3TBT3.fasta
wt_seq = "".join(["MPYSNLHPAIPRPRGHRSKYVALIFLVASLMILWVAKDPPNHTLKYLALHLASHELGLLL",
"KNLCCLAEELCHVQSRYQGSYWKAVRACLGCPIHCMAMILLSSYFYFLQNTADIYLSWMF",
"GLLVLYKSLSMLLGLQSLTPAEVSAVCEEKKLNVAHGLAWSYYIGYLRLILPGLQARIRM",
"FNQLHNNMLSGAGSRRLYILFPLDCGVPDNLSVVDPNIRFRDMLPQQNIDRAGIKNRVYS",
"NSVYEILENGQPAGVCILEYATPLQTLFAMSQDAKAGFSREDRLEQAKLFCRTLEEILED",
"VPESRNNCRLIVYQEPTDGNSFSLSQEVLRHIRQEEKEEVTMNAPMTSVAPPPSVLSQEP",
"RLLISGMDQPLPLRTDLI"])

# Generate all single mutants
mutants = generate_single_mutant_dict(wt_seq)

'''B6 predictions (H-2d series)'''
# Get all WT antigens
predicted_epitopes = predictor.predict_sequences(
   sequences={'WT':wt_seq},
   alleles={
       'B6': ["H2-K*b","H2-D*b"]
   },
   result="all",
   comparison_quantity="affinity",
   verbose=0)

# Do batchwise prediction
for batch in generate_groups(0, len(mutants), 5):
    print(batch)
        
    predictions = predictor.predict_sequences(
       sequences=dict(list(mutants.items())[i] for i in batch),
       alleles={
           'B6': ["H2-K*b","H2-D*b"]
       },
       result="all",
       comparison_quantity="affinity",
       verbose=0)
    
    # Add to data, remove duplicates from unmutated WT residues
    predicted_epitopes = pd.concat([predicted_epitopes,predictions]
                                   ).drop_duplicates(subset=['peptide',
                                                             'sample_name'],
                                                     ignore_index=True)

predicted_epitopes.to_csv('mhcflurry_outputs_B6.csv')

'''BALBc predictions (H-2b series)'''
# Get all WT antigens
predicted_epitopes = predictor.predict_sequences(
   sequences={'WT':wt_seq},
   alleles={
       "BALBc": ["H2-K*d","H2-L*d","H2-D*d"]
   },
   result="all",
   comparison_quantity="affinity",
   verbose=0)

# Do batchwise prediction
for batch in generate_groups(0, len(mutants), 5):
    print(batch)
        
    predictions = predictor.predict_sequences(
       sequences=dict(list(mutants.items())[i] for i in batch),
       alleles={
           "BALBc": ["H2-K*d","H2-L*d","H2-D*d"]
       },
       result="all",
       comparison_quantity="affinity",
       verbose=0)
    
    # Add to data, remove duplicates from unmutated WT residues
    predicted_epitopes = pd.concat([predicted_epitopes,predictions]
                                   ).drop_duplicates(subset=['peptide',
                                                             'sample_name'],
                                                     ignore_index=True)

predicted_epitopes.to_csv('mhcflurry_outputs_BALBc.csv')



