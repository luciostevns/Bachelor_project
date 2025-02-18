#%%

import pandas as pd
from Uniprot_ID_TO_organism_ID import with_organism_name
from Blosom_alignment import process_blosom_alignment
from allignment import process_alignment

if __name__ == "__main__":
    pathogen_DRB1_0401 = pd.read_csv("./Data/pathogen_DRB1_0401.csv")
    IEDB_DRB1_0401 = pd.read_csv("./Data/IEDB_DRB1_0401.csv")
    pathogen_DRB1_1501 = pd.read_csv("./Data/pathogen_DRB1_1501.csv")
    IEDB_DRB1_1501 = pd.read_csv("./Data/IEDB_DRB1_1501.csv")
    pathogen_DRB5_0101 = pd.read_csv("./Data/pathogen_DRB5_0101.csv")
    IEDB_DRB5_0101 = pd.read_csv("./Data/IEDB_DRB5_0101.csv")

    # Define pairs of DataFrames
    pairs = [
        (pathogen_DRB1_0401, IEDB_DRB1_0401),
        (pathogen_DRB1_1501, IEDB_DRB1_1501),
        (pathogen_DRB5_0101, IEDB_DRB5_0101)
    ]

    # Process the pairs
    combined_blosom = process_blosom_alignment(pairs, 2)
    combined_sim = process_alignment(pairs, 60, 90)

    blosom_finished = with_organism_name(combined_blosom, "./Data/uniprot_sprot.fasta")
    sim_finished = with_organism_name(combined_sim, "./Data/uniprot_sprot.fasta")

    # Save the combined DataFrame to a CSV file
    blosom_finished.to_csv("../Data/results/blosom_finished.csv", index=False)
    sim_finished.to_csv("../Data/results/sim_finished.csv", index=False)

# %%
