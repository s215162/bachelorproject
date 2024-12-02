import numpy as np
import pandas as pd
import logging

# Setup logging
log_file = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/matrix_weight_summation_GLCTLVAML"
)
logging.basicConfig(
    filename=log_file,
    filemode="a",  # Append mode
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO,
)
console = logging.StreamHandler()  # Log to console as well
console.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

# ESM file path
logging.info(f"Loading esm-embedding file")
esm_file = '/net/mimer/mnt/tank/projects2/emison/language_model/final_work/NEW_esm_embedding_matrix_full_data_GLCTLVAML.npz'

# Load the embeddings and raw indexes
esm_data = np.load(esm_file)
logging.info(f"Extracting embeddings and raw_indexes")
embeddings = esm_data['embeddings']  # Shape should be (num_rows, 1280)
raw_indexes = esm_data['raw_indexes']

# Load the CSV containing sequence data
logging.info(f"Loading csv info-file")
csv_file = '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv'
sequence_data = pd.read_csv(csv_file)

# Specify the peptide you want to process (e.g., 'GLCTLVAML')
peptide_x = 'GLCTLVAML'

# Filter the sequence data for the peptide
logging.info(f"Filter sequence data based on {peptide_x}")
peptide_data = sequence_data[sequence_data['peptide_x'] == peptide_x]

# List of CDR columns you need to sum
cdr_columns = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR1b', 'CDR2b', 'CDR3b']
logging.info(f"Working with weights {cdr_columns}")


# Loop over each CDR column (i.e., each "weight")
for cdr_column in cdr_columns:
    logging.info(f"Currently at weight {cdr_column}...")
    # Create a new list to hold summed embeddings for each raw_index
    summed_embeddings = []

    logging.info(f"Looping over each row to extract only the sequence of the weight")
    # Loop over each row in the peptide data
    for _, row in peptide_data.iterrows():
        raw_index = row['raw_index']  # Extract raw_index for this row
        if raw_index in raw_indexes:
            matrix_index = np.where(raw_indexes == raw_index)[0][0]  # Find corresponding matrix index
            matrix = embeddings[matrix_index]  # Get the embedding matrix for this raw_index

            # Extract the CDR positions for this CDR column
            cdr_positions = row[cdr_column].split(',')  # Assuming the positions are stored as comma-separated strings

            # Convert to integer positions
            cdr_positions = list(map(int, cdr_positions))
            logging.info(f"Summing embedding for weight {cdr_column}")
            # Sum the embeddings at these positions
            cdr_summed = np.zeros(matrix.shape[1])  # Initialize an array for the summed embeddings (1280 elements)
            for pos in cdr_positions:
                cdr_summed += matrix[pos]  # Add the corresponding row from the matrix

            summed_embeddings.append([raw_index] + list(cdr_summed))  # Store raw_index and summed embeddings

    logging.info(f"Convert data to numpy")
    # Convert the summed embeddings to a numpy array
    summed_embeddings = np.array(summed_embeddings)
    logging.info(f"Create dataframe for weight {cdr_column}")
    # Create a dataframe to save as CSV
    summed_embeddings_df = pd.DataFrame(summed_embeddings, columns=['raw_index'] + [f'embedding_{i}' for i in range(1280)])

    # Define the output CSV file name based on peptide_x and the current CDR column
    output_csv = f'/net/mimer/mnt/tank/projects2/emison/language_model/final_work/{peptide_x}_{cdr_column}_summed_embeddings.csv'
    logging.info(f"Saving data...")
    # Save the dataframe as a CSV file
    summed_embeddings_df.to_csv(output_csv, index=False)

    logging.info(f"Saved summed embeddings for {cdr_column} to {output_csv}.")