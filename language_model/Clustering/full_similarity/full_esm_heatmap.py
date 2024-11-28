import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

# # Directory containing the CSV files
# directory = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/redo'

# # List of files matching the pattern {peptide}_full.csv
# files = [f for f in os.listdir(directory) if f.endswith('_full.csv')]

# # List to hold dataframes
# dataframes = []

# # Iterate over the files
# for file in files:
#     # Read the current file
#     file_path = os.path.join(directory, file)
#     df = pd.read_csv(file_path)
    
#     # Extract peptide name from the file name (assuming the file is named {peptide}_full.csv)
#     peptide_name = file.split('_')[0]
    
#     # Add a new column 'peptide' with the peptide name
#     df['peptide'] = peptide_name
    
#     # Append the dataframe to the list
#     dataframes.append(df)

# # Concatenate all dataframes into one
# final_df = pd.concat(dataframes, ignore_index=True)

# Save the concatenated file with peptides included
output_file = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/concatenated_full_file.csv'
# final_df.to_csv(output_file, index=False)

# print(f"Concatenation complete! Saved to {output_file}")

# Load the concatenated data
data = pd.read_csv(output_file)

# Extract the 'sum_max_similarity' values for the pairwise similarity calculation
# You need to make sure you have 'raw_index', 'peptide', and 'sum_max_similarity' columns in your data
data_for_similarity = data.pivot_table(index='raw_index', columns='peptide', values='sum_max_similarity')

# Ensure that all columns contain numeric data (this is crucial for distance calculation)
data_for_similarity = data_for_similarity.apply(pd.to_numeric, errors='coerce')

# Handle any potential NaN values resulting from non-numeric columns, if necessary
# However, you mentioned no NaNs in your data, so this should be avoided.
data_for_similarity = data_for_similarity.fillna(0)

# Calculate the pairwise absolute differences
# Using broadcasting to subtract each pair of rows from one another
difference_matrix = data_for_similarity.values[:, None, :] - data_for_similarity.values[None, :, :]

# Calculate the absolute difference
difference_matrix = abs(difference_matrix)

# Convert the difference matrix to a DataFrame
difference_df = pd.DataFrame(difference_matrix.reshape(data_for_similarity.shape[0], -1), 
                             index=data_for_similarity.index, 
                             columns=data_for_similarity.columns)

# Save the difference matrix
output_file2 = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/difference_matrix.csv'
difference_df.to_csv(output_file2)

# Plot the heatmap of the difference matrix
plt.figure(figsize=(10, 8))
sns.heatmap(difference_df, annot=False, cmap='coolwarm', linewidths=0.5, square=True, cbar_kws={'label': 'Difference'})
plt.title('Pairwise Difference Matrix Heatmap')
plt.xlabel('Peptide')
plt.ylabel('Peptide')

# Save the heatmap
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/difference_matrix_heatmap.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save

print(f"Difference matrix saved to {output_file2}")
print(f"Heatmap saved to {output_path}")
