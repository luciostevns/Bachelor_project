#%%
import pandas as pd 
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Load the MMseqs2 clustering result
reps = pd.read_csv("../Data/clusterRes_cluster.tsv", sep="\t", header=None, names=["Representative", "Member"])

# Read the FASTA file
ref_prots_fasta = list(SeqIO.parse("../Data/clusterRes_rep_seq.fasta", "fasta"))
all_fasta = list(SeqIO.parse("../Data/all_proteomes.fasta", "fasta"))

def parse_fasta_to_df(fasta): 
    metadata = []
    
    for seq_record in fasta:
        header = seq_record.description
        
        # Extract the protein ID
        protein_id_match = re.search(r'\|([^|]+)\|', header)
        protein_id = protein_id_match.group(1) if protein_id_match else seq_record.id
        
        # Extract the organism source
        organism_source = header.split('OS=')[1].split(' OX=')[0] if "OS=" in header else None
        if organism_source:
            organism_source = ' '.join(organism_source.split()[:2])  # Keep only genus and species
        
        # Extract the taxonomy ID
        taxonomy_id_match = re.search(r'OX=(\d+)', header)
        taxonomy_id = taxonomy_id_match.group(1) if taxonomy_id_match else None
        
        # Extract the annotation (exclude the repeated protein ID)
        annotation = None
        header_parts = header.split()
        if len(header_parts) > 1:  # Ensure there is an annotation part
            possible_annotation = ' '.join(header_parts[1:])  # Everything after the first part
            annotation = possible_annotation.split('OS=')[0].strip()  # Stop at OS= if present
        
        # Add extracted data to list
        metadata.append([protein_id, organism_source, taxonomy_id, annotation, str(seq_record.seq)])
    
    # Convert to DataFrame
    metadata_df = pd.DataFrame(metadata, columns=["Protein ID", "Organism Source", "Taxonomy ID", "Annotation", "Sequence"])
    
    return metadata_df



# Run the function to retrieve protein metadata from fastas
rep_df = parse_fasta_to_df(ref_prots_fasta)
all_df = parse_fasta_to_df(all_fasta)

# Save to csv
all_df.to_csv("../Data/wrangled_all_pathogen_prots.csv", index=False)
rep_df.to_csv("../Data/wrangled_rep_pathogen_prots.csv", index=False)


#------------------------------ PLOTING -------------------------------#

# Define the datasets and labels
datasets = [("Representative Proteins", rep_df), ("All Proteins", all_df)]

# Define the columns to plot
columns = ["Organism Source"]

# Loop over datasets and columns
for dataset_name, df in datasets:
    for column in columns:
        if column in df.columns:  # Ensure the column exists
            counts = df[column].value_counts()

            # Plot bar chart (Top 20)
            plt.figure(figsize=(12, 6))
            sns.barplot(x=counts.index[:20], y=counts.values[:20], palette="viridis")
            plt.xticks(ticks=range(len(counts.index[:20])), labels=counts.index[:20], rotation=45, ha="right")
            plt.xlabel(column)
            plt.ylabel("Protein Count")
            plt.title(f"Distribution of {dataset_name} by {column} (Top 20)")
            plt.grid(axis="y", linestyle="--", alpha=0.7)
            plt.show()

# Count the number of members per representative protein (each cluster)
cluster_sizes = reps["Representative"].value_counts()

# Plot the distribution of cluster sizes
plt.figure(figsize=(10, 6))
plt.hist(cluster_sizes, bins=30, edgecolor="black", alpha=0.75)
plt.xlabel("Proteins per Cluster")
plt.ylabel("Number of Clusters")
plt.title("Distribution of Proteins per Cluster in MMseqs2 Output")
plt.yscale("log")  # Log scale helps visualize if distribution is highly skewed
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show plot
plt.show()
    
# %%
