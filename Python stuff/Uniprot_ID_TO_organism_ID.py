import pandas as pd

def with_organism_name(similarity_df, fasta_file_path):
    # Read the FASTA file
    with open(fasta_file_path, 'r') as file:
        fasta_content = file.read()

    # Split the content by '>' to separate each entry
    entries = fasta_content.split('>')[1:]

    # Initialize lists to store the UniProt IDs and organism sources
    uniprot_ids = []
    organism_sources = []

    # Iterate through each entry
    for entry in entries:
        # Extract the header line
        header_line = entry.split('\n')[0]
        
        # Extract the UniProt ID if it exists
        if '|' in header_line:
            uniprot_id = header_line.split('|')[1]
        else:
            uniprot_id = None
        
        # Extract the organism source if it exists
        if 'OS=' in header_line:
            organism_source = header_line.split('OS=')[1].split(' OX=')[0]
        else:
            organism_source = None
        
        # Append the extracted information to the lists
        uniprot_ids.append(uniprot_id)
        organism_sources.append(organism_source)

    # Create a DataFrame
    data = {'UniProt ID': uniprot_ids, 'Organism Source': organism_sources}
    df = pd.DataFrame(data)

    # Merge the DataFrames using a left join
    merged_df = similarity_df.merge(df, left_on='Pathogen_pep', right_on='UniProt ID', how='left')

    # Remove rows where Pathogen_pep is NA
    merged_df = merged_df.dropna(subset=['Pathogen_pep'])

    # Drop the UniProt ID column
    merged_df = merged_df.drop(columns=['UniProt ID'])

    return merged_df