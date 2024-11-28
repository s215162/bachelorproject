import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

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

# Pivot the data to create a matrix of sum_max_similarity for each peptide
# The index can be something like 'binder' or 'TCR ID', depending on the structure of your data
pivot_data = final_df.pivot_table(index='binder', columns='peptide', values='sum_max_similarity', aggfunc='mean')

# Optional: Normalize the data (if needed) using MinMaxScaler (scales values to [0, 1])
# scaler = MinMaxScaler()
# pivot_data_scaled = scaler.fit_transform(pivot_data)

# Convert back to a DataFrame for the heatmap
# pivot_data_scaled_df = pd.DataFrame(pivot_data_scaled, columns=pivot_data.columns, index=pivot_data.index)
pivot_data_df = pd.DataFrame(pivot_data, columns=pivot_data.columns, index=pivot_data.index)

# Generate a custom color palette for peptides
peptide_colors = sns.color_palette("Set2", len(pivot_data.columns))  # You can change the palette
peptide_color_map = dict(zip(pivot_data.columns, peptide_colors))

# Create a heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(
    pivot_data_df,
    annot=True,  # Annotate the cells with the data
    cmap=sns.color_palette(peptide_colors),  # Use the custom color palette
    fmt='.2f',  # Format values to two decimal places
    cbar_kws={'label': 'Similarity Score'},  # Add a color bar
    xticklabels=pivot_data.columns,  # Show peptide names on the x-axis
    yticklabels=pivot_data.index  # Show binder names on the y-axis
)

# Add a title and labels
plt.title('Heatmap of Peptide Data for sum of ESM2 values for full tcra+b sequence')
plt.xlabel('Peptide')
plt.ylabel('Binder')

# Save the heatmap as an image file
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/heatmap_peptide_clusters_unique_colors.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save
plt.show()
