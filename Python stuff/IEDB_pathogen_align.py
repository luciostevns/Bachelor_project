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
epitope_ids = IEDB_data["Assay_ID"].to_numpy()
epitope_sources = IEDB_data["Epitope - Molecule Parent"]

protein_sequences = pathogen_data["Sequence"].to_numpy()
protein_ids = pathogen_data["Protein ID"].to_numpy()
organism = pathogen_data["Organism Source"].to_numpy()

# Function to generate all 9-mers from an epitope sequence
def generate_9mers(seq):
    return [seq[i:i+9] for i in range(len(seq) - 8)] if len(seq) >= 9 else []

# Create a mapping of 9-mers to their respective Assay_IDs and Epitope Source
nine_mer_to_assay = {}
for assay_id, epitope_source, seq in zip(epitope_ids, epitope_sources, epitope_sequences):
    for nine_mer in generate_9mers(seq):
        # Now store both assay_id and epitope_source in the dictionary
        nine_mer_to_assay[nine_mer] = (assay_id, epitope_source)  # Store tuple of (Assay_ID, Epitope Source)

# Aho-Corasick Trie for fast multiple substring search
A = ahocorasick.Automaton()
for nine_mer, (assay_id, epitope_source) in nine_mer_to_assay.items():
    A.add_word(nine_mer, (assay_id, epitope_source, nine_mer))  # Store Assay_ID and Epitope Source with 9-mer
A.make_automaton()  # Build trie

# Store results
matches = []

# Function to find matches in a protein sequence
def find_matches(protein, protein_id, organism):
    local_matches = []
    for _, (assay_id, epitope_source, nine_mer) in A.iter(protein):
        local_matches.append((assay_id, epitope_source, protein_id, organism, nine_mer))
    return local_matches

# Iterate over protein sequences and store matches
with tqdm(total=len(protein_sequences), desc="Processing Proteins") as pbar:
    for protein, protein_id, organism in zip(protein_sequences, protein_ids, organism):
        matches.extend(find_matches(protein, protein_id, organism))
        pbar.update(1)

# Convert to DataFrame for easier viewing & saving
match_df = pd.DataFrame(matches, columns=["Assay_ID", "Epitope Source", "Protein_ID", "Organism Source", "Matched_9mer"])

# Print total matches and head of DataFrame
print(f"\nTotal Matches Found: {len(matches)}")


"""
# Plot distribution of matches per Assay_ID
plt.figure(figsize=(10, 6))
sns.countplot(y="Assay_ID", data=match_df, order=match_df["Assay_ID"].value_counts().index)
plt.title("Distribution of Matches per Assay_ID")
plt.xlabel("Number of Matches")
plt.ylabel("Assay_ID")
plt.tight_layout()
plt.show()

# Plot distribution of matches per Protein_ID
plt.figure(figsize=(10, 6))
sns.countplot(y="Protein_ID", data=match_df, order=match_df["Protein_ID"].value_counts().index)
plt.title("Distribution of Matches per Protein ID")
plt.xlabel("Number of Matches")
plt.ylabel("Protein ID")
plt.tight_layout()
plt.show()

# Plot distribution of matches per Organism Source
plt.figure(figsize=(10, 6))
sns.countplot(y="Organism Source", data=match_df, order=match_df["Organism Source"].value_counts().index)
plt.title("Distribution of Matches per Organism Source")
plt.xlabel("Number of Matches")
plt.ylabel("Organism Source")
plt.tight_layout()
plt.show()

"""
# %%
