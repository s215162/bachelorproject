import pandas as pd
import os
import sys

def merge_dataframes(df1, df2):
    """Merge two DataFrames on 'raw_index' and 'binder'."""
    return pd.merge(df1, df2, on=["raw_index", "binder", "partition"], how='outer')

def further_merge(peptide):
    """Merge two CSV datafiles with columns 'raw_index', 'sum_max_similarity', 'binder', and 'partition' and sum their max_similarity values."""
    # Define the directory where the files are located
    directory = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC"
    
    # Construct file paths for CDR1 and CDR2
    fp1 = os.path.join(directory, f"{peptide}_CDR1_CDR2_merged.csv")  # Adjust the file extension as needed
    fp2 = os.path.join(directory, f"{peptide}_CDR3.csv")  # Adjust the file extension as needed

    # Read the CDR files
    df1 = pd.read_csv(fp1)
    df2 = pd.read_csv(fp2)

    # Print columns for debugging
    print(f"File 1 Columns: {df1.columns.tolist()}")  # Debugging
    print(f"File 2 Columns: {df2.columns.tolist()}")  # Debugging

    # Merge CDR1 and CDR2 on 'raw_index' and 'binder'
    merged_cdr = merge_dataframes(df1, df2)

    # Sum the max similarities from CDR1 and CDR2
    merged_cdr['sum_max_similarity'] = (
        merged_cdr['sum_max_similarity_x'].fillna(0) + 
        merged_cdr['sum_max_similarity_y'].fillna(0)
    )

    # Create a final DataFrame with the desired columns
    final_cdr_df = merged_cdr[["raw_index", "sum_max_similarity", "binder", "partition"]]

    # Save to a new file
    final_filename = os.path.join(directory, f"{peptide}_Unweighed.csv")
    final_cdr_df.to_csv(final_filename, index=False)
    print(f"Merged file saved as {final_filename}")

# Check for command-line argument
if len(sys.argv) != 2:
    print("Usage: python unweighed.py <peptide>")
    sys.exit(1)

# Define your peptide from command-line argument
peptide = sys.argv[1]

# Execute the further merge
further_merge(peptide)
