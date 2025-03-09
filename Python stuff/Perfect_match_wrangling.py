import pandas as pd

# Load the data
perfect_match_rep = pd.read_csv("../Data/perfect_matches_rep.csv")
rep_pathogen_prot = pd.read_csv("../Data/wrangled_all_pathogen_prots.csv")

# Filter Protein_IDs where "Epitope Source" contains "Histone"
histone_proteins_ids = perfect_match_rep[perfect_match_rep["Epitope Source"].str.contains("Histone", na=False)]["Protein_ID"]

# Merge to get the sequences for these Protein_IDs
histone_proteins = rep_pathogen_prot[rep_pathogen_prot["Protein ID"].isin(histone_proteins_ids)]

# Write the sequences to a FASTA file for BLASTing
fasta_filename = "histone_proteins.fasta"
with open(fasta_filename, "w") as fasta_file:
    for _, row in histone_proteins.iterrows():
        protein_id = row["Protein ID"]
        sequence = row["Sequence"]
        fasta_file.write(f">{protein_id}\n{sequence}\n")

