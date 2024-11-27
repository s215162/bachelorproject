print("Importing...")
import numpy as np
import pandas as pd

print("Reading inputfile...")
# Read inputfile
# CSV file, each row represents a position in a TCR sequence. Consists of a raw_index column and the embedding column
input_file = '/net/mimer/mnt/tank/projects2/emison/language_model/ELAGIGILTV_tcr_matrices.csv'
output_file = '/net/mimer/mnt/tank/projects2/emison/language_model/ELAGIGILTV_tcr_VECTORS.csv'

print("Reading input data...")
# Read inputdata
data = pd.read_csv(input_file, low_memory=False)

print("group embeddings...")
# Sum each row
summed_embeddings = data.sum(axis=1)

print("Embeddings summed... saving file")
# Save as new CSV file
summed_embeddings.to_csv(output_file, index=True)  # 'index=True' to keep 'raw_index' in the output

print(f'Summed embeddings saved to {output_file}')
