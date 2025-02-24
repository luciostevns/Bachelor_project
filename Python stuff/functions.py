import pandas as pd 
from Bio.Align import substitution_matrices
from Bio import SeqIO
import re

def process_blosom_alignment(pairs, threshold):
    """
    Process BLOSUM alignment for pairs of DataFrames.

    Parameters:
    pairs (list of tuples): List of tuples where each tuple contains two DataFrames (pathogen_data, IEDB_data).
    threshold (float): Similarity threshold for comparing sequences.

    Returns:
    pd.DataFrame: Combined DataFrame containing similarity results.
    """        
    
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
        """
        Compare properties of two sequences using BLOSUM62 matrix.

        Parameters:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        threshold (float): Similarity threshold.

        Returns:
        float or None: Similarity score if above threshold, otherwise None.
        """

        similarity_results = []

        for _, p_row in pathogen_data.iterrows():
            for _, a_row in IEDB_data.iterrows():
                
                # Compare the properties
                core_sim = compare_properties(a_row["IEDB_core"], p_row["pathogen_core"], threshold)
                tcore_sim = compare_properties(a_row["IEDB_tcore"], p_row["pathogen_tcore"], threshold)
                
                if core_sim is not None and tcore_sim is not None:
                    combined_row = {**p_row.to_dict(), **a_row.to_dict()}
                    combined_row['core_sim'] = core_sim
                    combined_row['tcore_sim'] = tcore_sim
                    similarity_results.append(combined_row)

        similarity_df = pd.DataFrame(similarity_results)
        
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


def extract_cores_from_proteome(fasta_file_path, pattern_tcore_list, sample_size):
    """
    Extract cores from a proteome FASTA file based on a given pattern.

    Parameters:
    fasta_file_path (str): Path to the proteome FASTA file.
    patterns (list of str): Regular expression pattern to detect in the sequences.

    Returns:
    pd.DataFrames: DataFrames containing the Protein ID, Core sequences and organism source.
    """
    # Read the FASTA file
    proteome_fasta = list(SeqIO.parse(fasta_file_path, "fasta"))
    proteome_fasta = proteome_fasta[:sample_size]

    # Create a dictionary to store the results for each pattern
    results = {pattern: [] for pattern, _ in pattern_tcore_list}

    # Iterate over the DataFrame and search for the pattern in each sequence
    for seq_record in proteome_fasta:
        header = seq_record.description
        sequence = str(seq_record.seq)
        
        for pattern, tcore_list in pattern_tcore_list:
            # Find all matches of the pattern in the sequence
            matches = re.finditer(pattern, sequence)
        
            for match in matches:
                pathogen_core = match.group()
                pathogen_tcore = ''.join([pathogen_core[i-1] for i in tcore_list])
                protein_id_match = re.search(r'\|([^|]+)\|', header)
                organism_source = header.split('OS=')[1].split(' OX=')[0] if "OS=" in header else None
                if protein_id_match:
                    protein_id = protein_id_match.group(1)
                    results[pattern].append([protein_id, pathogen_core, pathogen_tcore, organism_source])

    # Convert the results to a DataFrame
    cores_dfs = {pattern: pd.DataFrame(data, columns=["Protein ID", "pathogen_core", "pathogen_tcore", "Organism_source"]) for pattern, data in results.items()}

    return cores_dfs

