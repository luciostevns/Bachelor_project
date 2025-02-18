import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Align import substitution_matrices

def process_blosom_alignment(pairs, threshold):
    blosum62 = substitution_matrices.load("BLOSUM62")

    def compare_properties(seq1, seq2, threshold):
        total_similarity = 0
        for a1, a2 in zip(seq1, seq2):
            total_similarity += blosum62[a1][a2]
        
        if total_similarity / len(seq1) > threshold:
            return total_similarity / len(seq1)
        else:
            return None

    def get_similarity(pathogen_data, IEDB_data):
        similarity_results = []

        for _, p_row in pathogen_data.iterrows():
            for _, a_row in IEDB_data.iterrows():
                
                # Compare the properties
                core_sim = compare_properties(a_row["core"], p_row["core"], threshold)
                tcore_sim = compare_properties(a_row["core"], p_row["core"], threshold)
                
                if core_sim is not None and tcore_sim is not None:
                    similarity_results.append([p_row['Pathogen_prot'], p_row["Seq"], a_row['Epitope_ID'], a_row["Seq"], core_sim, tcore_sim])

        similarity_df = pd.DataFrame(similarity_results, columns=['Pathogen_pep', "Pathogen_seq", 'Epitope_ID', "IEDB_seq", 'core_sim', "tcore_sim"])
        
        return similarity_df

    combined_results = []

    for pathogen_df, IEDB_df in pairs:
        # Get similarity results
        blosom_results = get_similarity(pathogen_df, IEDB_df)

        # Remove duplicates
        unique_blosom_results = blosom_results.drop_duplicates()

        combined_results.append(unique_blosom_results)

    # Concatenate all results
    combined_df = pd.concat(combined_results, axis=0)

    return combined_df