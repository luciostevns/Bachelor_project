#%%
import numpy as np
from Bio import pairwise2
from Bio.Align import substitution_matrices
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

# Define binary hydropathy (0 = hydrophobic, 1 = hydrophilic)
amino_acid_properties = {
    'A': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'R': {'hydropathy': 1, 'charge': +1, 'aromatic': 0},  # Hydrophilic
    'N': {'hydropathy': 1, 'charge': 0, 'aromatic': 0},  # Hydrophilic
    'D': {'hydropathy': 1, 'charge': -1, 'aromatic': 0},  # Hydrophilic
    'C': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'Q': {'hydropathy': 1, 'charge': 0, 'aromatic': 0},  # Hydrophilic
    'E': {'hydropathy': 1, 'charge': -1, 'aromatic': 0},  # Hydrophilic
    'G': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'H': {'hydropathy': 1, 'charge': +1, 'aromatic': 0},  # Hydrophilic
    'I': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'L': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'K': {'hydropathy': 1, 'charge': +1, 'aromatic': 0},  # Hydrophilic
    'M': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'F': {'hydropathy': 0, 'charge': 0, 'aromatic': 1},  # Hydrophobic
    'P': {'hydropathy': 0, 'charge': 0, 'aromatic': 0},  # Hydrophobic
    'S': {'hydropathy': 1, 'charge': 0, 'aromatic': 0},  # Hydrophilic
    'T': {'hydropathy': 1, 'charge': 0, 'aromatic': 0},  # Hydrophilic
    'W': {'hydropathy': 1, 'charge': 0, 'aromatic': 1},  # Hydrophilic
    'Y': {'hydropathy': 1, 'charge': 0, 'aromatic': 1},  # Hydrophilic
    'V': {'hydropathy': 0, 'charge': 0, 'aromatic': 0}   # Hydrophobic
}


# Read the CSV file into a DataFrame
pathogen_DRB1_0401 = pd.read_csv("pathogen_DRB1_0401.csv")
IEDB_DRB1_0401 = pd.read_csv("IEDB_DRB1_0401.csv")
pathogen_DRB1_1501 = pd.read_csv("pathogen_DRB1_1501.csv")
IEDB_DRB1_1501 = pd.read_csv("IEDB_DRB1_1501.csv")
pathogen_DRB5_0101 = pd.read_csv("pathogen_DRB5_0101.csv")
IEDB_DRB5_0101 = pd.read_csv("IEDB_DRB5_0101.csv")

def sequence_to_properties(sequence):
    return [amino_acid_properties[aa] for aa in sequence]

def compare_properties(seq1_properties, seq2_properties, threshold):
    total_similarity = 0
    for p1, p2 in zip(seq1_properties, seq2_properties):
        # Compare binary hydropathy, charge, and aromatic for each position
        hydropathy_sim = 1 if p1['hydropathy'] == p2['hydropathy'] else 0
        charge_sim = 1 if p1['charge'] == p2['charge'] else 0
        aromatic_sim = 1 if p1['aromatic'] == p2['aromatic'] else 0
        
        # Sum the similarities
        similarity = hydropathy_sim + charge_sim + aromatic_sim
        total_similarity += similarity
    
    similarity_percentage = (total_similarity / (len(seq1_properties) * 3)) * 100
    if similarity_percentage > threshold:
        return similarity_percentage
    else:
        return None


def get_similarity(pathogen_data, IEDB_data):
    similarity_results = []

    for _, p_row in pathogen_data.iterrows():
        for _, a_row in IEDB_data.iterrows():
            # Convert core sequences to properties
            p_core_prop = sequence_to_properties(p_row['core'])
            a_core_prop = sequence_to_properties(a_row['core'])

            tp_core_prop = sequence_to_properties(p_row['TCR_bind_core'])
            ta_core_prop = sequence_to_properties(a_row['TCR_bind_core'])
            
            # Compare the properties
            core_sim = compare_properties(p_core_prop, a_core_prop, 60)
            tcore_sim = compare_properties(tp_core_prop, ta_core_prop, 90)
            
            if core_sim is not None and tcore_sim is not None:
                similarity_results.append([p_row['Pathogen_prot'], p_row["Seq"], a_row['Epitope_ID'], a_row["Seq"], core_sim, tcore_sim])

    similarity_df = pd.DataFrame(similarity_results, columns=['Pathogen_pep', "Pathogen_seq", 'Epitope_ID', "IEDB_seq", 'core_sim', "tcore_sim"])
    
    return similarity_df


DRB1_0401_sim = get_similarity(pathogen_DRB1_0401, IEDB_DRB1_0401)
DRB1_1501_sim = get_similarity(pathogen_DRB1_1501, IEDB_DRB1_1501)
DRB5_0101_sim = get_similarity(pathogen_DRB5_0101, IEDB_DRB5_0101)

# Show the result
DRB1_0401_sim.to_csv("C:/Users/lucio/Desktop/DTU stuff/Bachelor stuff/Bachelor R stuff/Bachelor_project/Data/similarity_results_DRB1_0401.csv", index=False)
DRB1_1501_sim.to_csv("C:/Users/lucio/Desktop/DTU stuff/Bachelor stuff/Bachelor R stuff/Bachelor_project/Data/similarity_results_DRB1_1501.csv", index=False)
DRB5_0101_sim.to_csv("C:/Users/lucio/Desktop/DTU stuff/Bachelor stuff/Bachelor R stuff/Bachelor_project/Data/similarity_results_DRB5_0101.csv", index=False)

# Plot histograms
plt.figure(figsize=(12, 6))

# Histogram for core_sim
plt.subplot(1, 2, 1)
sns.histplot(DRB1_0401_sim['core_sim'], kde=True, color='blue', bins=30)
plt.title('Histogram of Core Similarity (DRB1_0401)')
plt.xlabel('Core Similarity (%)')
plt.ylabel('Frequency')

# Histogram for tcore_sim
plt.subplot(1, 2, 2)
sns.histplot(DRB1_0401_sim['tcore_sim'], kde=True, color='green', bins=30)
plt.title('Histogram of TCR Core Similarity (DRB1_0401)')
plt.xlabel('TCR Core Similarity (%)')
plt.ylabel('Frequency')

plt.tight_layout()
plt.show()

# %%
