import torch
import esm
import pandas as pd
import os
import logging
import sys
import numpy as np  # Added for saving .npz files

# Command-line input for the peptide
if len(sys.argv) != 2:
    print("Usage: python script.py <peptide>")
    sys.exit(1)

peptide = sys.argv[1]

# Setup logging
log_file = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/TEST_full_ESM_"
    + peptide
    + "_matrix_cdr_processing.log"
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

# Define the dataset path
dataset_path = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
)

# Define output directory
output_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # Disables dropout for deterministic results

# Load and filter dataset based on the specified peptide
logging.info(f"Loading and filtering dataset for peptide {peptide}")
df = pd.read_csv(dataset_path)
df_filtered = df[df["peptide_x"] == peptide]


# Function to process and save sequence representations
def process_sequences(df_filtered, peptide):
    sequence_representations = []
    raw_indexes = []
    batch_size = 5  # Adjust based on memory capacity
    for start in range(0, batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        data = [(row["peptide_x"], row["tcr_full"]) for _, row in df_subset.iterrows()]
        raw_indexes.extend(df_subset["raw_index"].tolist())  # Collect raw_index values
        logging.info(f"Raw-indexes are: {raw_indexes}")
        logging.info(
            f"Processing rows {start} to {end} for peptide {peptide}, column {'tcr_full'}"
        )

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            # Append the full-length 2D embeddings
            sequence_representations.append(token_representations[i, 1 : tokens_len - 1].cpu().numpy())

    logging.info(f"ESM finished.. Saving file...")

    # Save as .npz to preserve full 2D embeddings
    npz_file_path = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/TEST_esm_embedding_matrix_full_data_{peptide}.npz"
    np.savez(npz_file_path, embeddings=sequence_representations, raw_indexes=raw_indexes)
    logging.info(f"Saved 2D embeddings to {npz_file_path}")

# Process and save summed vectors for the specified dataset as
logging.info(f"Processing dataset for peptide {peptide}")
process_sequences(df_filtered, peptide)
