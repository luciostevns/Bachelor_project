import pandas as pd

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

def sequence_to_properties(sequence):
    return [amino_acid_properties[aa] for aa in sequence]

def process_alignment(pairs, core_threshold, tcore_threshold):
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
                core_sim = compare_properties(p_core_prop, a_core_prop, core_threshold)
                tcore_sim = compare_properties(tp_core_prop, ta_core_prop, tcore_threshold)
                
                if core_sim is not None and tcore_sim is not None:
                    similarity_results.append([p_row['Pathogen_prot'], p_row["Seq"], a_row['Epitope_ID'], a_row["Seq"], core_sim, tcore_sim])

        similarity_df = pd.DataFrame(similarity_results, columns=['Pathogen_pep', "Pathogen_seq", 'Epitope_ID', "IEDB_seq", 'core_sim', "tcore_sim"])
        
        return similarity_df

    combined_results = []

    for pathogen_df, IEDB_df in pairs:
        # Get similarity results
        alignment_results = get_similarity(pathogen_df, IEDB_df)

        # Remove duplicates
        unique_alignment_results = alignment_results.drop_duplicates()

        combined_results.append(unique_alignment_results)

    # Concatenate all results
    combined_df = pd.concat(combined_results, axis=0)

    return combined_df