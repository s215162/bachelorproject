print("Importing...")
import numpy as np
import pandas as pd

peptides = [
        "RAKFQLL",
        "YLQPRTFLL",
]

for peptide in peptides:

    print(f"Reading input file for peptide {peptide}...")
    input_file = f'/net/mimer/mnt/tank/projects2/emison/language_model/binders_matrix/{peptide}_tcr_matrices.csv'
    output_file = f'/net/mimer/mnt/tank/projects2/emison/language_model/binders_matrix/{peptide}_tcr_VECTORS.csv'

    print("Reading input data...")
    data = pd.read_csv(input_file, low_memory=False)

    # Check the shape of the data
    print(f"Data shape: {data.shape}")
    
    # Check the first few rows to ensure it's loading correctly
    print("Summing embeddings...")

    # Sum each row (excluding non-embedding columns)
    embedding_cols = data.columns[2:]  # 'raw_index' and 'sequence' are the first two columns
    summed_embeddings = data[embedding_cols].sum(axis=1)

    # Print the first few values (for debugging)
    print(summed_embeddings[:5])

    print("Embeddings summed... saving file")
    
    # Save as new CSV file (including 'raw_index')
    output_data = pd.DataFrame({
        'raw_index': data['raw_index'],
        'summed_embeddings': summed_embeddings
    })

    output_data.to_csv(output_file, index=False)

    print(f'Summed embeddings saved to {output_file}')

