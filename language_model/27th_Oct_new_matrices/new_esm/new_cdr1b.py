import torch
import esm
import pandas as pd
import os
import logging
import sys

# Command-line input for the peptide
if len(sys.argv) != 2:
    print("Usage: python script.py <peptide>")
    sys.exit(1)

peptide = sys.argv[1]

# Setup logging
log_file = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/logs/cdr1b_esm_processing.log"
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

# Define the datasets
datasets = {
    "binders": "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv",
    "swaps": "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv",
}

# Define output directory for CDR1 beta
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/beta/cdr1"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Filter datasets based on the specified peptide to reduce data size upfront
logging.info(f"Loading and filtering datasets for peptide {peptide}")
df_binders = pd.read_csv(datasets["binders"])
df_filtered_binders = df_binders[df_binders["peptide_x"] == peptide]

df_swaps = pd.read_csv(datasets["swaps"])
df_filtered_swaps = df_swaps[df_swaps["peptide_x"] == peptide]

# Function to process and save sequence representations
def process_sequences(df_filtered, peptide, dataset_name):
    all_sequence_representations = []
    raw_indexes = []
    batch_size = 25  # Adjust based on memory capacity
    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        data = [(row["peptide_x"], row["CDR1b"]) for _, row in df_subset.iterrows()]
        raw_indexes.extend(df_subset["raw_index"].tolist())  # Collect raw_index values
        logging.info(f"Processing rows {start} to {end} for {dataset_name}, peptide {peptide}")

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            embedding_matrix = token_representations[i, 1 : tokens_len - 1].cpu().numpy()
            all_sequence_representations.append(embedding_matrix)

    # Create DataFrame with raw_index and flattened embeddings for export
    output_df = pd.DataFrame(
        {"raw_index": raw_indexes, "embedding_matrix": all_sequence_representations}
    )

    output_file_path = os.path.join(
        output_dir, f"sequence_embeddings_{dataset_name}_{peptide}.csv"
    )
    output_df.to_csv(output_file_path, index=False)
    logging.info(f"Saved to {output_file_path}")

# Process and save embeddings for binders and swaps
logging.info(f"Processing binders for peptide {peptide}")
process_sequences(df_filtered_binders, peptide, "binders")

logging.info(f"Processing swaps for peptide {peptide}")
process_sequences(df_filtered_swaps, peptide, "swaps")

