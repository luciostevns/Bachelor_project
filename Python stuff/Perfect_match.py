#%%

import pandas as pd
import numpy as np
import ahocorasick
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
pathogen_data = pd.read_csv("./Data/wrangled_rep_pathogen_prots.csv")
IEDB_data = pd.read_csv("./Data/wrangled_IEDB.csv")

# Extract relevant columns
epitope_sequences = IEDB_data["Sequence"].to_numpy()
assay_id = IEDB_data["Assay_ID"].to_numpy()
epitope_sources = IEDB_data["Protein_source"].to_numpy()
disease = IEDB_data["Disease"].to_numpy()

protein_sequences = pathogen_data["Sequence"].to_numpy()
protein_ids = pathogen_data["Protein ID"].to_numpy()
organism = pathogen_data["Organism Source"].to_numpy()
annotation = pathogen_data["Annotation"].to_numpy()
location = pathogen_data["Subcellular location [CC]"].to_numpy()

# Function to generate all 9-mers from an epitope sequence
def generate_9mers(seq):
    return [seq[i:i+9] for i in range(len(seq) - 8)] if len(seq) >= 9 else []

# Create a mapping of 9-mers to their respective Assay_IDs and Epitope Source
nine_mer_to_assay = {}
for assay_id, epitope_source, disease, seq in zip(assay_id, epitope_sources, disease, epitope_sequences):
    for nine_mer in generate_9mers(seq):
        # Now store both assay_id and epitope_source in the dictionary
        nine_mer_to_assay[nine_mer] = (assay_id, epitope_source, disease)  # Store tuple of (Assay_ID, Epitope Source)

# Aho-Corasick Trie for fast multiple substring search
A = ahocorasick.Automaton()
for nine_mer, (assay_id, epitope_source, disease) in nine_mer_to_assay.items():
    A.add_word(nine_mer, (assay_id, epitope_source, disease, nine_mer))  # Store Assay_ID and Epitope Source with 9-mer
A.make_automaton()  # Build trie

# Store results
matches = []

# Function to find matches in a protein sequence
def find_matches(protein, protein_id, organism, annotation, location):
    local_matches = []
    for _, (assay_id, epitope_source, disease, nine_mer) in A.iter(protein):
        local_matches.append((assay_id, epitope_source, disease, protein_id, organism, annotation, location, nine_mer))
    return local_matches

# Iterate over protein sequences and store matches
with tqdm(total=len(protein_sequences), desc="Processing Proteins") as pbar:
    for protein, protein_id, organism, annotation, location in zip(protein_sequences, protein_ids, organism, annotation, location):
        matches.extend(find_matches(protein, protein_id, organism, annotation, location))
        pbar.update(1)

# Convert to DataFrame for easier viewing & saving
match_df = pd.DataFrame(matches, columns=["Assay_ID", "Epitope Source", "Disease", "Protein_ID", "Organism Source", "Pathogen Annotation", "Subcellular location (pathogen)", "Matched_9mer"])

match_df.to_csv("./Data/perfect_matches_rep.csv", index=False)

# Print total matches and head of DataFrame
print(f"\nTotal Matches Found: {len(matches)}")

#--------------------------- PLOTING ------------------------------#

columns_to_plot = [
    "Assay_ID",
    "Epitope Source",
    "Disease",
    "Protein_ID",
    "Organism Source",
    "Pathogen Annotation"
]

for col in columns_to_plot:
    plt.figure(figsize=(10, 6))
    sns.countplot(y=col, data=match_df, order=match_df[col].value_counts().head(20).index)
    plt.title(f"Distribution of Matches per {col}")
    plt.xlabel("Number of Matches")
    plt.ylabel(col)
    plt.tight_layout()
    plt.show()

# %%
