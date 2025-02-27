#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
    
# Pathogen reference data
pathogen_ref = pd.read_csv("./Data/pathogen_ref.tsv", sep='\t')

# Proteome data
referenced_proteome = pd.read_csv("./Data/referenced_proteomes.tsv", sep='\t')
other_proteome = pd.read_csv("./Data/other_proteomes.tsv", sep='\t')

# binding rows of ref and other
proteome = pd.concat([referenced_proteome, other_proteome])

# Define the pattern you want to detect (e.g., 'human' or 'homo' in any case)
pattern = r'(?i)human|homo'

# Remove instances where the 'Host' column is NaN
pathogen_ref = pathogen_ref.dropna(subset=['Host'])

pathogen_ref = pathogen_ref.dropna(subset=['Assembly'])

# Filter out rows that don't contain the pattern in the 'Host' column
filtered_df = pathogen_ref[pathogen_ref['Host'].str.contains(pattern, regex=True)]

# Exclude rows that contain 'non-human' in the 'Host' column
filtered_pathogen_ref = filtered_df[~filtered_df['Host'].str.contains(r'(?i)non', regex=True)]

# Merging proteome data onto ref
# Look into why more rows is occuring, should only have 1 ID
merged_proteome = filtered_pathogen_ref.merge(proteome, how="left", left_on="Assembly", right_on="Genome assembly ID")

# Remove weird instances with missing Proteime id
# Looks at NA's
unique_merged = merged_proteome.dropna(subset=['Proteome Id'])

# Saving certain columns from merged df as .csv
unique_merged[["Proteome Id", "#Organism group","Strain", "Protein count", "Assembly"]].to_csv("./Data/Pathogenic_bacteria_proteome.csv", index=False)

# Saving proteomes IDs as .txt for download of full proteomes
unique_merged['Proteome Id'].to_csv("./Data/proteome_ids.txt", header=False, index=False)

print(unique_merged["Protein count"].sum())

# %%
