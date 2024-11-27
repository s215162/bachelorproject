import pandas as pd

# Load the TCR data file that contains the 'binder' column and other details
tcr_info_df = pd.read_csv(
    "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
)

peptide = "ELAGIGILTV"

# Filter only rows where peptide matches "ELAGIGILTV"
tcr_info_df = tcr_info_df[tcr_info_df["peptide_x"] == peptide]

# Define peptide and file paths for alpha chain files
to_add_binder_alpha = [
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/full/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_full_alpha.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr1/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr1a.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr2/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr2a.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr3/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr3a.csv",
]

count = 0

for data in to_add_binder_alpha:
    # Load the max_similarity data
    similarity_df = pd.read_csv(data)

    # Merge the dataframes on the 'raw_index' column, ensuring only ELAGIGILTV peptide rows are included
    merged_df = pd.merge(
        similarity_df, tcr_info_df[["raw_index", "binder"]], on="raw_index", how="left"
    )

    # Determine the file label based on count
    if count == 0:
        dis = "full"
    elif count == 1:
        dis = "cdr1"
    elif count == 2:
        dis = "cdr2"
    else:
        dis = "cdr3"

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(
        f"data_with_binder_" + peptide + "_alpha_" + dis + ".csv", index=False
    )

    count += 1

# Repeat the same steps for beta chain files
to_add_binder_beta = [
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/full/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_full_beta.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr1/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr1b.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr2/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr2b.csv",
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr3/"
    + peptide
    + "/max_cosine_similarity_results_"
    + peptide
    + "_cdr3b.csv",
]

count = 0

for data in to_add_binder_beta:
    # Load the max_similarity data
    similarity_df = pd.read_csv(data)

    # Merge the dataframes on the 'raw_index' column, ensuring only ELAGIGILTV peptide rows are included
    merged_df = pd.merge(
        similarity_df, tcr_info_df[["raw_index", "binder"]], on="raw_index", how="left"
    )

    # Determine the file label based on count
    if count == 0:
        dis = "full"
    elif count == 1:
        dis = "cdr1"
    elif count == 2:
        dis = "cdr2"
    else:
        dis = "cdr3"

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(
        f"data_with_binder_" + peptide + "_beta_" + dis + ".csv", index=False
    )

    count += 1

