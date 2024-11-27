import torch
import esm
import pandas as pd
import os
import logging
import sys

# Command-line input for the peptide and weight
if len(sys.argv) != 3:
    print("Usage: python script.py <peptide> <weight>")
    sys.exit(1)

peptide = sys.argv[1]
weight = sys.argv[2]  # This will be the weight type, e.g., CDR2a

# Setup logging
log_file = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/logs/cdr_processing_binders.log"
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

# Define the dataset path for binders only
dataset_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv"

# Define output directory based on weight type
output_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/{weight.lower()}"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Load and filter dataset based on the specified peptide to reduce data size upfront
logging.info(f"Loading and filtering binders dataset for peptide {peptide}")
df = pd.read_csv(dataset_path)
df_filtered = df[df["peptide_x"] == peptide]

# Function to process and save sequence representations
def process_sequences(df_filtered, peptide, weight, dataset_name):
    all_sequence_representations = []
    raw_indexes = []
    batch_size = 25  # Adjust based on memory capacity
    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        # Use the weight variable to select the corresponding column
        data = [(row["peptide_x"], row[weight]) for _, row in df_subset.iterrows()]
        raw_indexes.extend(df_subset["raw_index"].tolist())  # Collect raw_index values
        logging.info(f"Processing rows {start} to {end} for {dataset_name}, peptide {peptide}, weight {weight}")

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            embedding_matrix = token_representations[i, 1 : tokens_len - 1].cpu().numpy()
            all_sequence_representations.append(embedding_matrix)

    # Create DataFrame with raw_index and embeddings for export
    output_df = pd.DataFrame(
        {"raw_index": raw_indexes, "embedding_matrix": all_sequence_representations}
    )

    output_file_path = os.path.join(
        output_dir, f"sequence_embeddings_binders_{peptide}.csv"
    )
    output_df.to_csv(output_file_path, index=False)
    logging.info(f"Saved to {output_file_path}")

# Process and save embeddings for the binders dataset
logging.info(f"Processing binders dataset for peptide {peptide} with weight {weight}")
process_sequences(df_filtered, peptide, weight, "binders")

