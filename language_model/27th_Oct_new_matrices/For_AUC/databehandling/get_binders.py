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
to_add_binder_alpha = [
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/full/{peptide}/max_cosine_similarity_results_{peptide}_full_alpha.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr1/{peptide}/max_cosine_similarity_results_{peptide}_cdr1a.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr2/{peptide}/max_cosine_similarity_results_{peptide}_cdr2a.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr3/{peptide}/max_cosine_similarity_results_{peptide}_cdr3a.csv",
]

count = 0

for data in to_add_binder_alpha:
    # Check for existence of the primary file
    if not os.path.exists(data):
        # Check for alternate filename with peptide in place of full_alpha
        alternative_data = data.replace("_full_alpha.csv", f"_full_{peptide}_alpha.csv")
        if os.path.exists(alternative_data):
            data = alternative_data
        else:
            # Check for the format that retains directory structure
            alternative_data = data.replace(".csv", f"_{peptide}.csv")
            if os.path.exists(alternative_data):
                data = alternative_data
            else:
                # Check for the format that retains directory structure
                alternative_data = data.replace("_alpha.csv", f"_{peptide}.csv")
                if os.path.exists(alternative_data):
                    data = alternative_data
                else:
                    print(f"File not found: {data} or {alternative_data}")
                    continue

    # Load the max_similarity data
    similarity_df = pd.read_csv(data)

    # Merge the dataframes on the 'raw_index' column
    merged_df = pd.merge(
        similarity_df, tcr_info_df[["raw_index", "binder"]], on="raw_index", how="left"
    )

    # Determine the file label based on count
    dis = ["full", "cdr1", "cdr2", "cdr3"][count]

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(f"data_with_binder_{peptide}_alpha_{dis}.csv", index=False)

    count += 1

# Repeat the same steps for beta chain files
to_add_binder_beta = [
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/full/{peptide}/max_cosine_similarity_results_{peptide}_full_beta.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr1/{peptide}/max_cosine_similarity_results_{peptide}_cdr1b.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr2/{peptide}/max_cosine_similarity_results_{peptide}_cdr2b.csv",
    f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr3/{peptide}/max_cosine_similarity_results_{peptide}_cdr3b.csv",
]

count = 0

for data in to_add_binder_beta:
    # Check for existence of the primary file
    if not os.path.exists(data):
        # Check for alternate filename with peptide in place of full_beta
        alternative_data = data.replace("_beta.csv", f"_{peptide}_beta.csv")
        if os.path.exists(alternative_data):
            data = alternative_data
        else:
            # Check for the format that retains directory structure
            alternative_data = data.replace(".csv", f"_{peptide}.csv")
            if os.path.exists(alternative_data):
                data = alternative_data
            else:
                # Check for the format that retains directory structure
                alternative_data = data.replace("_beta.csv", f"_{peptide}.csv")
                if os.path.exists(alternative_data):
                    data = alternative_data
                else:
                    # Check for the format that retains directory structure
                    alternative_data = data.replace("_beta.csv", f"_{peptide}.csv")
                    if os.path.exists(alternative_data):
                        data = alternative_data
                    else:
                            print(f"File not found: {data} or {alternative_data}")
                            continue

    # Load the max_similarity data
    similarity_df = pd.read_csv(data)

    # Merge the dataframes on the 'raw_index' column
    merged_df = pd.merge(
        similarity_df, tcr_info_df[["raw_index", "binder"]], on="raw_index", how="left"
    )

    # Determine the file label based on count
    dis = ["full", "cdr1", "cdr2", "cdr3"][count]

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(f"data_with_binder_{peptide}_beta_{dis}.csv", index=False)

    count += 1

