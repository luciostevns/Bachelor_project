#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read the TSV files
pathogen_ref = pd.read_csv("./Data/pathogen_ref.tsv", sep='\t')
proteome = pd.read_csv("./Data/proteome.tsv", sep='\t')

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
merged_proteome = filtered_pathogen_ref.merge(proteome, how="left", left_on="Assembly", right_on="Genome assembly ID")

# Remove weird instances with missing Proteime id
unique_merged = merged_proteome.dropna(subset=['Proteome Id'])

# Saving certain columns from merged df as .csv
unique_merged[["Proteome Id", "#Organism group","Strain", "Protein count", "Assembly"]].to_csv("./Data/Pathogenic_bacteria_proteome.csv", index=False)

print("unique: ", unique_merged.shape)
print("merged: ", merged_proteome.shape)
# %%
