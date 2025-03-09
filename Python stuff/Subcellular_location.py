#%%
import pandas as pd

rep_df = pd.read_csv("./Data/wrangled_rep_pathogen_prots.csv")
all_df = pd.read_csv("./Data/wrangled_all_pathogen_prots.csv")
cell_location_ref = pd.read_csv("./Data/cell_location_uniprot.tsv", sep="\t")

# Function to extract the first location
cell_location_ref["Subcellular location [CC]"] = cell_location_ref["Subcellular location [CC]"].str.extract(r"SUBCELLULAR LOCATION:\s*([^;{]+)")

# Assigning cellular location to proteins
merged_rep_df = rep_df.merge(cell_location_ref, how="left", left_on="Protein ID", right_on="Entry")
merged_all_df = all_df.merge(cell_location_ref, how="left", left_on="Protein ID", right_on="Entry")

# Save to CSV
#merged_rep_df.to_csv("./Data/wrangled_rep_pathogen_prots.csv", index=False)
#merged_all_df.to_csv("./Data/wrangled_all_pathogen_prots.csv", index=False)

na_percent_rep = merged_rep_df["Subcellular location [CC]"].isna().mean() * 100
na_percent_all = merged_all_df["Subcellular location [CC]"].isna().mean() * 100

print(f"Percentage of missing values in merged_rep_df: {na_percent_rep:.2f}%")
print(f"Percentage of missing values in merged_all_df: {na_percent_all:.2f}%")
# %%
