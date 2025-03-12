#%%

import pandas as pd
import numpy as np
import ahocorasick
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
rep_pathogen_data = pd.read_csv("../Data/wrangled_rep_pathogen_prots.csv")
all_pathogen_data = pd.read_csv("../Data/wrangled_all_pathogen_prots.csv")

IEDB_data = pd.read_csv("../Data/wrangled_IEDB.csv")

def find_epitope_matches(pathogen_data, IEDB_data):
    # Extract relevant columns, including 'Subcellular location [CC]'
    epitope_sequences = IEDB_data["Sequence"].to_numpy()
    assay_id = IEDB_data["Assay_ID"].to_numpy()
    epitope_sources = IEDB_data["Protein_source"].to_numpy()
    disease = IEDB_data["Disease"].to_numpy()

    protein_sequences = pathogen_data["Sequence"].to_numpy()
    protein_ids = pathogen_data["Protein ID"].to_numpy()
    organism = pathogen_data["Organism Source"].to_numpy()
    annotation = pathogen_data["Annotation"].to_numpy()
    subcellular_location = pathogen_data["Subcellular location [CC]"].to_numpy()  # New column

    # Function to generate all 9-mers from an epitope sequence
    def generate_9mers(seq):
        return [seq[i:i+9] for i in range(len(seq) - 8)] if len(seq) >= 9 else []

    # Create a mapping of 9-mers to their respective Assay_IDs and Epitope Source
    nine_mer_to_assay = {}
    for assay_id, epitope_source, disease, seq in zip(assay_id, epitope_sources, disease, epitope_sequences):
        for nine_mer in generate_9mers(seq):
            nine_mer_to_assay[nine_mer] = (assay_id, epitope_source, disease)

    # Aho-Corasick Trie for fast multiple substring search
    A = ahocorasick.Automaton()
    for nine_mer, (assay_id, epitope_source, disease) in nine_mer_to_assay.items():
        A.add_word(nine_mer, (assay_id, epitope_source, disease, nine_mer))
    A.make_automaton()

    # Store results
    matches = []

    # Function to find matches in a protein sequence
    def find_matches(protein, protein_id, organism, annotation, subcellular_location):
        local_matches = []
        for _, (assay_id, epitope_source, disease, nine_mer) in A.iter(protein):
            local_matches.append((assay_id, epitope_source, disease, protein_id, organism, annotation, subcellular_location, nine_mer))
        return local_matches

    # Iterate over protein sequences and store matches
    with tqdm(total=len(protein_sequences), desc="Processing Proteins") as pbar:
        for protein, protein_id, organism, annotation, subcellular_location in zip(protein_sequences, protein_ids, organism, annotation, subcellular_location):
            matches.extend(find_matches(protein, protein_id, organism, annotation, subcellular_location))
            pbar.update(1)

    # Convert to DataFrame for easier viewing & saving
    match_df = pd.DataFrame(matches, columns=["Assay_ID", "Epitope Source", "Disease", "Protein_ID", "Organism Source", "Pathogen Annotation", "Subcellular Location", "Matched_9mer"])
    
    return match_df


rep_df = find_epitope_matches(rep_pathogen_data, IEDB_data)
all_df = find_epitope_matches(all_pathogen_data, IEDB_data)

rep_df.to_csv("../Data/perfect_matches_rep.csv", index=False)
all_df.to_csv("../Data/perfect_matches_all.csv", index=False)

# --------------------------- PLOTTING ------------------------------#

def plot_match_distribution(match_df, dataset_name):
    """Plots the distribution of matches for selected columns."""
    columns_to_plot = [
        "Assay_ID",
        "Epitope Source",
        "Disease",
        "Protein_ID",
        "Organism Source",
        "Pathogen Annotation",
        "Subcellular Location"  # Added this column to the plot
    ]

    for col in columns_to_plot:
        plt.figure(figsize=(10, 6))

        # Check if the column has at least one non-null value
        if match_df[col].notna().sum() == 0:
            print(f"Skipping plot for {col} as it contains only NaN values.")
            continue

        # Plot the top 20 most frequent values
        sns.countplot(y=col, data=match_df, order=match_df[col].value_counts().head(20).index)
        plt.title(f"Distribution of Matches per {col} ({dataset_name})")
        plt.xlabel("Number of Matches")
        plt.ylabel(col)
        plt.tight_layout()
        plt.show()

# Call the function for each dataset
plot_match_distribution(rep_df, "Rep Pathogen Data")
plot_match_distribution(all_df, "All Pathogen Data")

print("Percentage of matches in Rep Pathogen Data:", len(rep_df) / len(rep_pathogen_data) * 100)
print("Number of rep proteins with matches:", len(rep_df["Protein_ID"]))
print("Number of unique rep proteins with matches:", len(rep_df["Protein_ID"].unique()))

print("Percentage of matches in All Pathogen Data:", len(all_df) / len(all_pathogen_data) * 100)
print("Number of all proteins with matches:", len(all_df["Protein_ID"]))
print("Number of unique all proteins with matches:", len(all_df["Protein_ID"].unique()))

# %%
