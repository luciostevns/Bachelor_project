#%%

import pandas as pd

# Load the data
perfect_match_rep = pd.read_csv("../Data/perfect_matches_rep.csv")
perfect_match_all = pd.read_csv("../Data/perfect_matches_all.csv")
rep_pathogen_prot = pd.read_csv("../Data/wrangled_all_pathogen_prots.csv")
all_pathogen_prot = pd.read_csv("../Data/wrangled_all_pathogen_prots.csv")


def filter_and_export_fasta(match_df, pathogen_df, pattern, output_fasta):
    # Filter Protein_IDs where "Epitope Source" contains the given pattern
    filtered_protein_ids = match_df[match_df["Pathogen Annotation"].str.contains(pattern, na=False)]["Protein_ID"]
    
    # Retrieve sequences for the filtered Protein_IDs
    filtered_proteins = pathogen_df[pathogen_df["Protein ID"].isin(filtered_protein_ids)]
    
    # Write sequences to a FASTA file
    with open(output_fasta, "w") as fasta_file:
        for _, row in filtered_proteins.iterrows(): 
            fasta_file.write(f">{row['Protein ID']}\n{row['Sequence']}\n")

filter_and_export_fasta(perfect_match_rep, rep_pathogen_prot, "Histone", "../Data/histone_proteins.fasta")
filter_and_export_fasta(perfect_match_all, all_pathogen_prot, "Chaperone protein DnaK", "../Data/Chap_DnaK_proteins.fasta")


# %%
