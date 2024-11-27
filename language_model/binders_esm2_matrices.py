import pandas as pd
import torch
import esm
import os
import logging

# Configure logging to save to a file
log_file_path = '/net/mimer/mnt/tank/projects2/emison/language_model/binders_matrix/process.log'
logging.basicConfig(
    filename=log_file_path,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger()

# Logging the start of the process
logger.info("Importing modules...")
print("Importing...")

logger.info("Loading model")
print("Loading model")
# Load the ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # Disable dropout for deterministic results

peptide_names = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFQLL",
    "YLQPRTFLL",
]

logger.info("Defining file paths")
print("Stating file paths")
# File paths
input_file = '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv'
output_dir = '/net/mimer/mnt/tank/projects2/emison/language_model/binders_matrix'

logger.info("File paths defined")
print("File paths defined")

# Load the dataset
logger.info("Loading data from CSV")
df = pd.read_csv(input_file)
logger.info("Data loaded")
print("Data loaded")

for peptide in peptide_names:
    output_file = os.path.join(output_dir, peptide+'_tcr_matrices.csv')
    # Filter for rows where 'peptide_x' is 'peptide'
    filtered_df = df[df['peptide_x'] == peptide]
    logger.info(f"Data filtered for peptide: {peptide}")
    print("Data filtered")

    # List to store results for each sequence
    results_list = []

    # Process each sequence in 'tcr_full'
    for index, row in filtered_df.iterrows():
        logger.info(f"Processing index {index}")
        print("Processing index {}".format(index))
        tcr_full_seq = row['tcr_full']

        # Prepare data for the model
        data = [("tcr_full", tcr_full_seq)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)

        logger.info("Data prepared for model")
        print("Data prepared")

        # Extract per-residue representations (on CPU)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)
        token_representations = results["representations"][33]
        logger.info("Per-residue representations extracted")
        print("per-residue representations extracted")

        # Remove special tokens (beginning and end) and convert to a matrix
        seq_len = (batch_tokens != alphabet.padding_idx).sum(1).item()
        matrix = token_representations[0, 1:seq_len-1].numpy()  # Extract the matrix as a numpy array

        logger.info("Converted to matrix")
        print("Converted to matrix")

        # Append to the results list
        results_list.append({
            "raw_index": row["raw_index"],
            "sequence": tcr_full_seq,
            "matrix": matrix
        })
        logger.info("Appended results for index {}".format(index))
        print("Appended to results")

    # Save results to CSV (matrix as string representation)
    logger.info("Saving results to CSV")
    print("Saving file")
    with open(output_file, 'w') as f_out:
        f_out.write("raw_index,sequence,matrix\n")
        for result in results_list:
            matrix_str = '\n'.join(','.join(map(str, row)) for row in result['matrix'])
            f_out.write(f"{result['raw_index']},{result['sequence']},{matrix_str}\n")

    logger.info(f"Results saved to {output_file}")
    print(f"Results saved to {output_file}")

