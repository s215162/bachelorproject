import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import pairwise_distances

# Directory containing the CSV files
directory = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/redo'

# List of files matching the pattern {peptide}_full.csv
files = [f for f in os.listdir(directory) if f.endswith('_full.csv')]

# List to hold dataframes
dataframes = []

# Iterate over the files
for file in files:
    # Read the current file
    file_path = os.path.join(directory, file)
    df = pd.read_csv(file_path)
    
    # Extract peptide name from the file name (assuming the file is named {peptide}_full.csv)
    peptide_name = file.split('_')[0]
    
    # Add a new column 'peptide' with the peptide name
    df['peptide'] = peptide_name
    
    # Append the dataframe to the list
    dataframes.append(df)

# Concatenate all dataframes into one
final_df = pd.concat(dataframes, ignore_index=True)

# Save the concatenated file with peptides included
output_file = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/concatenated_full_file.csv'
final_df.to_csv(output_file, index=False)

print(f"Concatenation complete! Saved to {output_file}")

# Load the concatenated data
data = pd.read_csv(output_file)

# Extract the 'sum_max_similarity' values for the pairwise similarity calculation
# You need to make sure you have 'raw_index', 'peptide', and 'sum_max_similarity' columns in your data
data_for_similarity = data.pivot_table(index='raw_index', columns='peptide', values='sum_max_similarity')

# Ensure that we have the peptide columns correctly formatted
# You may want to normalize your data (if desired) using MinMaxScaler
scaler = MinMaxScaler()
data_for_similarity_scaled = scaler.fit_transform(data_for_similarity)

# Calculate the pairwise similarity (or distance) between rows
similarity_matrix = pairwise_distances(data_for_similarity_scaled, metric='euclidean')  # Use 'cosine' or 'manhattan' if needed

# Convert the similarity matrix back to a DataFrame
similarity_df = pd.DataFrame(similarity_matrix, index=data_for_similarity.index, columns=data_for_similarity.index)
output_file2 = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/sim_matrix.csv'
final_df.to_csv(output_file2, index=False)

# Plot the heatmap of the similarity matrix
plt.figure(figsize=(10, 8))
sns.heatmap(similarity_df, annot=False, cmap='coolwarm', linewidths=0.5, square=True, cbar_kws={'label': 'Similarity (scaled)'})
plt.title('Pairwise Similarity Matrix Heatmap')
plt.xlabel('Raw Index')
plt.ylabel('Raw Index')

# Save the heatmap
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/similarity_matrix_heatmap.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save