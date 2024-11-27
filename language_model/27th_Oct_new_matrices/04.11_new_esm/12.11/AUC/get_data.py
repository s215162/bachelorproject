import pandas as pd
import os
import sys

# Check for command-line argument
if len(sys.argv) != 2:
    print("Usage: python get_data.py <peptide>")
    sys.exit(1)

peptide = sys.argv[1]

def read_and_rename(file, suffix):
    """Read a CSV file and rename the 'max_similarity' column."""
    df = pd.read_csv(file)
    return df.rename(columns={"max_similarity": f"max_similarity_{suffix}"})

def merge_dataframes(df1, df2):
    """Merge two DataFrames on 'raw_index', 'binder', and 'partition'."""
    return pd.merge(df1, df2, on=["raw_index", "binder", "partition"])

def calculate_sum_max_similarity(merged_df):
    """Calculate the sum of max similarities."""
    return merged_df["max_similarity_file1"] + merged_df["max_similarity_file2"]

def process_files(peptide):
    """Process the files based on the weights."""
    weights = ["CDR1", "CDR2", "CDR3", "full"]
    filetypes_a = ["cdr1a", "cdr2a", "cdr3a", "tcra"]
    filetypes_b = ["cdr1b", "cdr2b", "cdr3b", "tcrb"]

    for weight in weights:
        alpha_type = weights.index(weight)
        beta_type = alpha_type

        # Load files
        file1 = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/similarities_with_binder_partition_{peptide}_{filetypes_a[alpha_type]}.csv"
        file2 = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/similarities_with_binder_partition_{peptide}_{filetypes_b[beta_type]}.csv"

        # Read and merge
        df1 = read_and_rename(file1, "file1")
        df2 = read_and_rename(file2, "file2")
        merged_df = merge_dataframes(df1, df2)

        # Calculate sum max similarity
        merged_df["sum_max_similarity"] = calculate_sum_max_similarity(merged_df)

        # Save the result
        final_df = merged_df[["raw_index", "sum_max_similarity", "binder", "partition"]]
        output_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/{peptide}_{weight}.csv"
        final_df.to_csv(output_file, index=False)

process_files(peptide)
