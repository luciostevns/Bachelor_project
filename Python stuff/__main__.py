#%%

import pandas as pd
from functions import process_blosom_alignment
from functions import extract_cores_from_proteome
import time

if __name__ == "__main__":

    # Read in IEDB allele specific data
    IEDB_DRB1_0401 = pd.read_csv("./Data/IEDB_DRB1_0401.csv")
    IEDB_DRB1_1501 = pd.read_csv("./Data/IEDB_DRB1_1501.csv")
    IEDB_DRB5_0101 = pd.read_csv("./Data/IEDB_DRB5_0101.csv")

    # Define the patterns, ORDER IMPORTANT for later workflow insert
    # patterns in order of alleles
    pattern_tcore_list = [
        (r'[FYILV].{4}[ASGP].{2}[KR]', [2,3,4,5,7,8]),
        (r'[ILV].{2}[YF].{1}[SNGDA].{2}[LVAI]', [2,3,5,7,8])
    ]
    
    start_time = time.time()

    # Read in pathogen allele specific data
    pathogen_wrangled = extract_cores_from_proteome("./Data/all_proteomes.fasta", pattern_tcore_list, 100)
    
    print("time after pathogen wrangeling: ", time.time() - start_time)

    # Assign the DataFrames to variables with meaningful names
    pathogen_DRB1_0401 = pathogen_wrangled[pattern_tcore_list[0][0]]
    pathogen_DRB1_1501 = pathogen_wrangled[pattern_tcore_list[1][0]]

    # Define pairs of DataFrames
    pairs = [
        (pathogen_DRB1_0401, IEDB_DRB1_0401),
        (pathogen_DRB1_1501, IEDB_DRB1_1501)
    ]

    # Process the pairs
    combined_blosom = process_blosom_alignment(pairs, 3)

    print(f"Time aftr scorring: {time.time() - start_time} seconds")

    # Save the combined DataFrame to a CSV file
    combined_blosom.to_csv("../Data/results/blosom_finished.csv", index=False)


# %%
