import pandas as pd

gmt_file = "mouse_hallmarks_cGAMP_extended.gmt"
output_csv = "mouse_hallmarks_cGAMP_extended.csv"

data = []
with open(gmt_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        pathway_name = parts[0]
        description = parts[1]
        genes = parts[2:]
        for gene in genes:
            data.append([pathway_name, description, gene])

# Create DataFrame and export
df = pd.DataFrame(data, columns=["Pathway", "Description", "Gene"])
df.to_csv(output_csv, index=False)
