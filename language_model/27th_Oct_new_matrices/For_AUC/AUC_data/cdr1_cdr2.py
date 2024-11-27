import pandas as pd
import os

def merge_dataframes(df1, df2):
    """Merge two DataFrames on 'raw_index' and 'binder'."""
    return pd.merge(df1, df2, on=["raw_index", "binder"], how='outer')

def further_merge(peptide):
    """Merge CDR1 and CDR2 and sum their max_similarity values."""
    # Define the directory where the files are located
    directory = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/AUC_data"
    
    # Construct file paths for CDR1 and CDR2
    cdr1_file = os.path.join(directory, f"{peptide}_CDR1.csv")  # Adjust the file extension as needed
    cdr2_file = os.path.join(directory, f"{peptide}_CDR2.csv")  # Adjust the file extension as needed

    # Read the CDR files
    cdr1_df = pd.read_csv(cdr1_file)
    cdr2_df = pd.read_csv(cdr2_file)

    # Print columns for debugging
    print(f"CDR1 Columns: {cdr1_df.columns.tolist()}")  # Debugging
    print(f"CDR2 Columns: {cdr2_df.columns.tolist()}")  # Debugging

    # Merge CDR1 and CDR2 on 'raw_index' and 'binder'
    merged_cdr = merge_dataframes(cdr1_df, cdr2_df)

    # Sum the max similarities from CDR1 and CDR2
    merged_cdr['sum_max_similarity'] = (
        merged_cdr['sum_max_similarity_x'].fillna(0) + 
        merged_cdr['sum_max_similarity_y'].fillna(0)
    )

    # Create a final DataFrame with the desired columns
    final_cdr_df = merged_cdr[["raw_index", "sum_max_similarity", "binder"]]

    # Save to a new file
    final_filename = os.path.join(directory, f"{peptide}_CDR1_CDR2_merged.csv")
    final_cdr_df.to_csv(final_filename, index=False)
    print(f"Merged file saved as {final_filename}")

# Define your peptide
peptide = "ELAGIGILTV"

# Execute the further merge
further_merge(peptide)

