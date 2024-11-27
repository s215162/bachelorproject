import pandas as pd
import os

def read_and_rename(file, suffix):
    """Read a CSV file and rename the 'max_similarity' column."""
    df = pd.read_csv(file)
    print(f"Columns in {file} before renaming: {df.columns.tolist()}")  # Debugging
    return df.rename(columns={"max_similarity": f"max_similarity_{suffix}"})

def merge_dataframes(df1, df2):
    """Merge two DataFrames on 'raw_index' and 'binder'."""
    return pd.merge(df1, df2, on=["raw_index", "binder"])

def calculate_sum_max_similarity(merged_df, weight):
    """Calculate the sum of max similarities based on the weight."""
    if weight == "Weighed":
        return (
            merged_df["max_similarity_file1"]
            + merged_df["max_similarity_file2"]
            + 4 * merged_df["max_similarity_file3"]
        )
    else:
        return (
            merged_df["max_similarity_file1"]
            + merged_df["max_similarity_file2"]
            + merged_df["max_similarity_file3"]
        )

def process_files(peptide, weight):
    """Process the files based on the weight."""
    if weight in ["full", "CDR1", "CDR2", "CDR3"]:
        alpha_type = weights.index(weight)
        beta_type = alpha_type
        
        # Load files
        file1 = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/{peptide}/data_with_binder_{peptide}_{filetypes_a[alpha_type]}.csv"
        file2 = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/{peptide}/data_with_binder_{peptide}_{filetypes_b[beta_type]}.csv"
        outputname = f"{peptide}_{weight}.csv"

        # Read and merge
        df1 = read_and_rename(file1, "file1")
        df2 = read_and_rename(file2, "file2")
        merged_df = merge_dataframes(df1, df2)

        # Print merged DataFrame columns
        print(f"Merged columns (full/CDR): {merged_df.columns.tolist()}")  # Debugging

        # Calculate sum max similarity
        merged_df["sum_max_similarity"] = merged_df["max_similarity_file1"] + merged_df["max_similarity_file2"]

        # Save the result
        final_df = merged_df[["raw_index", "sum_max_similarity", "binder"]]
        final_df.to_csv(outputname, index=False)
        print(f"File saved as {outputname}")

    elif weight in ["Weighed", "Unweighed"]:
        # Load CDR files
        files = [
            f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/AUC_data/{peptide}_CDR1.csv",
            f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/AUC_data/{peptide}_CDR2.csv",
            f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/AUC_data/{peptide}_CDR3.csv"
        ]

        # Read and merge the first two CDR files
        df1 = read_and_rename(files[0], "file1")
        df2 = read_and_rename(files[1], "file2")
        merged_df = merge_dataframes(df1, df2)

        # Print merged DataFrame columns
        print(f"Merged columns (CDR part 1): {merged_df.columns.tolist()}")  # Debugging

        # Save the merged DataFrame as a temporary file
        temp_filename = f"temp_{peptide}.csv"
        merged_df.to_csv(temp_filename, index=False)

        # Read the temporary file and merge with the third CDR file
        df3 = read_and_rename(files[2], "file3")
        final_merged_df = merge_dataframes(pd.read_csv(temp_filename), df3)

        # Print merged DataFrame columns before calculation
        print(f"Merged columns (final): {final_merged_df.columns.tolist()}")  # Debugging

        # Calculate the summed max similarity
        final_merged_df["sum_max_similarity"] = calculate_sum_max_similarity(final_merged_df, weight)

        # Save the result
        final_df = final_merged_df[["raw_index", "sum_max_similarity", "binder"]]
        final_df.to_csv(f"{peptide}_{weight}.csv", index=False)
        print(f"File saved as {peptide}_{weight}.csv")

        # Optionally, remove the temporary file
        os.remove(temp_filename)

# Define your parameters
peptide = "ELAGIGILTV"
weights = ["full", "CDR1", "CDR2", "CDR3", "Unweighed", "Weighed"]
filetypes_a = ["alpha_full", "alpha_cdr1", "alpha_cdr2", "alpha_cdr3"]
filetypes_b = ["beta_full", "beta_cdr1", "beta_cdr2", "beta_cdr3"]

# Process each weight
for weight in weights:
    process_files(peptide, weight)

