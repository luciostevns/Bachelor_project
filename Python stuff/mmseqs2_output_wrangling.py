#%%
import pandas as pd 
from Bio import SeqIO
import re

# Load the MMseqs2 clustering result
df = pd.read_csv("./Data/clusterRes_cluster.tsv", sep="\t", header=None, names=["Representative", "Member"])

# Read the FASTA file
ref_prots_fasta = list(SeqIO.parse("./Data/clusterRes_rep_seq.fasta", "fasta"))

rep_organism = []

# Iterate over the sampled sequences and search for the patterns
for seq_record in ref_prots_fasta:
    header = seq_record.description
    
    # Extract the protein ID and Organism source
    protein_id_match = re.search(r'\|([^|]+)\|', header)
    protein_id = protein_id_match.group(1) if protein_id_match else seq_record.id
    organism_source = header.split('OS=')[1].split(' OX=')[0] if "OS=" in header else None
    
    # Shorten organism name to genus and species only
    if organism_source:
        organism_source = ' '.join(organism_source.split()[:2])  # Keep only genus and species

    # Add it to the list
    rep_organism.append([protein_id, organism_source, str(seq_record.seq)])

# Convert to df
rep_organism_df = pd.DataFrame(rep_organism, columns=["Protein ID", "Organism Source", "Sequence"])

# Save to CSV
rep_organism_df.to_csv("./Data/wrangled_rep_pathogen_prots.csv", index=False)



"""
plots:

# Count occurrences of each organism
organism_counts = df_protein_organism["Organism Source"].value_counts()

# Plot bar chart
plt.figure(figsize=(12, 6))
sns.barplot(x=organism_counts.index[:20], y=organism_counts.values[:20], palette="viridis")  # Show top 20 organisms

# Remove x-axis labels
plt.xticks([], [])  # Empty labels

# Labels and title
plt.xlabel("Organism Source")
plt.ylabel("Protein Count")
plt.title("Distribution of Proteins by Organism Source (Top 20)")
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show the plot
plt.show()

# Count the number of members per representative protein (each cluster)
cluster_sizes = df["Representative"].value_counts()

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
"""
    
# %%
