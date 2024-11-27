import pandas as pd

# Load your data from the CSV files
file1_path = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/ELAGIGILTV/data_with_binder_ELAGIGILTV_alpha_full.csv'  
file2_path = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/ELAGIGILTV/data_with_binder_ELAGIGILTV_beta_full.csv'  


file1 = pd.read_csv(file1_path)
file2 = pd.read_csv(file2_path)

# Print columns to confirm structure
print("File1 columns:", file1.columns)
print("File2 columns:", file2.columns)

# Merge the dataframes on 'raw_index'
merged_df = pd.merge(file1, file2, on="raw_index", suffixes=('_file1', '_file2'))

# Print columns after merge
print("Columns after merge:", merged_df.columns)

# Create max similarity column
similarity_columns_file1 = [col for col in merged_df.columns if 'similarity' in col and '_file1' in col]
if similarity_columns_file1:  # Ensure there are similarity columns
    merged_df["max_similarity_file1"] = merged_df[similarity_columns_file1].max(axis=1)
else:
    print("No similarity columns found in file1.")

# Check the new column's creation
print("Data after max similarity calculation:")
print(merged_df.head())

# Save the merged DataFrame to a new CSV file
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/debugging_file.csv' 
merged_df.to_csv(output_path, index=False)

print(f"Merged data saved to: {output_path}")

