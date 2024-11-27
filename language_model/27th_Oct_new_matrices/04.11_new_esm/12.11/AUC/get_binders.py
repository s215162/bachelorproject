import pandas as pd
import os
import sys

# Check for command-line argument
if len(sys.argv) != 2:
    print("Usage: python get_binders.py <peptide>")
    sys.exit(1)

peptide = sys.argv[1]

# Load the TCR data file that contains the 'binder' column and other details
tcr_info_df = pd.read_csv(
    "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
)

# Filter only rows where peptide matches the specified peptide
tcr_info_df = tcr_info_df[tcr_info_df["peptide_x"] == peptide]

# Define peptide and file paths for alpha chain files
to_add_binder_data = [
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr1a_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr1b_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr2a_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr2b_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr3a_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/cdr3b_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/tcra_{peptide}_similarity_results.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/tcrb_{peptide}_similarity_results.csv",

]

count = 0

for data in to_add_binder_data:
    # Check for existence of the primary file
    if not os.path.exists(data):
        print(f"File not found: {data}")
        continue

    # Load the max_similarity data
    similarity_df = pd.read_csv(data)

    # Merge the dataframes on the 'raw_index' column
    merged_df = pd.merge(
        similarity_df, tcr_info_df[["raw_index", "binder", "partition"]], on="raw_index", how="left"
    )

    # Determine the file label based on count
    dis = ["cdr1a", "cdr1b", "cdr2a", "cdr2b", "cdr3a", "cdr3b", "tcra", "tcrb"][count]

    # Save the merged dataframe to a new CSV file
    output_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/check/similarities_with_binder_partition_{peptide}_{dis}.csv"
    merged_df.to_csv(output_file, index=False)

    count += 1
