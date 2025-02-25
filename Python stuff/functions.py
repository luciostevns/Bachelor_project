import pandas as pd
from Bio.Align import substitution_matrices
from Bio import SeqIO
import re
import numpy as np
from itertools import product
import sys
import random

def process_blosom_alignment(pairs, threshold):
    """
    Process BLOSUM alignment for pairs of DataFrames without parallel processing.
    
    Parameters:
    pairs (list of tuples): List of (pathogen_data, IEDB_data) DataFrame pairs.
    threshold (float): Similarity threshold for filtering results.

    Returns:
    pd.DataFrame: Combined DataFrame containing similarity results.
    """
    
    ## **Step 1: Precompute BLOSUM62 Matrix for Fast Lookups**
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Create amino acid index mapping for fast lookup
    aa_to_index = {aa: i for i, aa in enumerate("ARNDCQEGHILKMFPSTWYVBZX*")}
    
    # Convert BLOSUM62 to a NumPy matrix
    blosum_matrix = np.zeros((26, 26), dtype=np.int32)
    for aa1, aa2 in product(aa_to_index.keys(), repeat=2):
        blosum_matrix[aa_to_index[aa1], aa_to_index[aa2]] = blosum62[aa1][aa2]
    
    def sequence_similarity(seq1, seq2):
        """ Compute sequence similarity using precomputed BLOSUM62 matrix. """
        indices1 = np.array([aa_to_index[aa] for aa in seq1])
        indices2 = np.array([aa_to_index[aa] for aa in seq2])
        similarity_score = np.sum(blosum_matrix[indices1, indices2])  
        return similarity_score / len(seq1)

    def compare_pair(p_row, a_row, threshold):
        """ Compare one pair of pathogen and IEDB sequences. """
        core_sim = sequence_similarity(a_row["IEDB_core"], p_row["pathogen_core"])
        tcore_sim = sequence_similarity(a_row["IEDB_tcore"], p_row["pathogen_tcore"])
        
        if core_sim > threshold and tcore_sim > threshold:
            return {**p_row.to_dict(), **a_row.to_dict(), 'core_sim': core_sim, 'tcore_sim': tcore_sim}
        return None

    def get_similarity_sequential(pathogen_data, IEDB_data, threshold):
        """ Compare sequences sequentially without parallel processing. """
        results = []

        # Iterate over each combination of rows
        for _, p_row in pathogen_data.iterrows():
            for _, a_row in IEDB_data.iterrows():
                result = compare_pair(p_row, a_row, threshold)
                if result:
                    results.append(result)

        return pd.DataFrame(results)
    
    # Process All Pairs and Merge Results
    combined_results = []
    for pathogen_df, IEDB_df in pairs:
        # Get similarity results using sequential processing
        blosum_results = get_similarity_sequential(pathogen_df, IEDB_df, threshold)
        
        # Remove duplicates early
        unique_blosum_results = blosum_results.drop_duplicates()
        
        combined_results.append(unique_blosum_results)

    return pd.concat(combined_results, axis=0) if combined_results else pd.DataFrame()


def extract_cores_from_proteome(fasta_file_path, pattern_tcore_list, sample_size):
    """
    Extract cores from a proteome FASTA file based on given patterns.

    Parameters:
    fasta_file_path (str): Path to the proteome FASTA file.
    pattern_tcore_list (list of tuples): List of tuples where each tuple contains a pattern (str) and tcore positions (list of int).
    sample_size (int): Number of sequences to sample from the FASTA file.

    Returns:
    dict of pd.DataFrame: Dictionary where keys are patterns and values are DataFrames containing the Protein ID, Core sequences, TCore sequences, and Organism source.
    """
    # Read the FASTA file
    proteome_fasta = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Sample the sequences
    sampled_sequences = random.sample(proteome_fasta, min(sample_size, len(proteome_fasta)))

    # Create a dictionary to store the results for each pattern
    results = {pattern: [] for pattern, _ in pattern_tcore_list}

    # Iterate over the sampled sequences and search for the patterns
    for seq_record in sampled_sequences:
        header = seq_record.description
        sequence = str(seq_record.seq)
        
        for pattern, tcore_positions in pattern_tcore_list:
            # Find all matches of the pattern in the sequence
            matches = re.finditer(pattern, sequence)
        
            for match in matches:
                pathogen_core = match.group()
                pathogen_tcore = ''.join([pathogen_core[i] for i in tcore_positions if i < len(pathogen_core)])
                protein_id_match = re.search(r'\|([^|]+)\|', header)
                organism_source = header.split('OS=')[1].split(' OX=')[0] if "OS=" in header else None
                if protein_id_match:
                    protein_id = protein_id_match.group(1)
                    results[pattern].append([protein_id, pathogen_core, pathogen_tcore, organism_source])

    # Convert the results to DataFrames
    cores_dfs = {pattern: pd.DataFrame(data, columns=["Protein ID", "pathogen_core", "pathogen_tcore", "Organism_source"]) for pattern, data in results.items()}

    return cores_dfs