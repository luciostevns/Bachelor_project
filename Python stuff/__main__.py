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
        (r'[ILV].{2}[YFW].{1}[SNGDA].{3}', [2,3,5,7,8]),
        (r'[YFIW].{2}[DEAIQ].{1}[TNSDV].{2}[ASGQ]', [2,3,5,7,8]),
        (r'[FYILVW].{4}[ASGPN].{2}[KR]', [2,3,4,5,7,8]),
    ]

    start_time = time.time()

    # Read in pathogen allele specific data
    pathogen_wrangled = extract_cores_from_proteome("./Data/all_proteomes.fasta", pattern_tcore_list, 2000)

    print(f"Time after core extraction: {time.time() - start_time} seconds")

    # Assign the DataFrames to variables with meaningful names
    pathogen_DRB1_1501 = pathogen_wrangled[pattern_tcore_list[0][0]]
    pathogen_DRB1_0401 = pathogen_wrangled[pattern_tcore_list[1][0]]
    pathogen_DRB5_0101 = pathogen_wrangled[pattern_tcore_list[2][0]]

    # Define pairs of DataFrames
    pairs = [
        (pathogen_DRB1_0401, IEDB_DRB1_0401),
        (pathogen_DRB1_1501, IEDB_DRB1_1501),
        (pathogen_DRB5_0101, IEDB_DRB5_0101)
    ]

    start_time = time.time()

    combined_blosom = process_blosom_alignment(pairs, 0)

    print(f"Time after blosom: {time.time() - start_time} seconds")

    # Save the combined DataFrame to a CSV file
    combined_blosom.to_csv("../Data/results/blosom_finished.csv", index=False)
# %%