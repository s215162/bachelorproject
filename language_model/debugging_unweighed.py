import torch
import esm
import pandas as pd
import os
import logging

# Set up logging
log_file = "/net/mimer/mnt/tank/projects2/emison/language_model/processing_log.txt"
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,  # You can use DEBUG for more verbosity
    format='%(asctime)s - %(levelname)s - %(message)s',
)

# Additionally log to console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# Define peptides and partitions
peptides = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFKQLL",
    "YLQPRTFLL",
]
partitions = [0, 1, 2, 3, 4]

# Define the datasets
datasets = {
    "binders_only": "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv",
    "full_sequence": "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv",
}

# Define output directory
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/esm_unweighed"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
logging.info("Loading ESM-2 model...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results
logging.info("Model loaded successfully.")


# Function to process and save sequence representations
def process_sequences(df_filtered, peptide, partition_desc, dataset_name):
    all_sequence_representations = []
    batch_size = 25  # Adjust based on memory capacity

    if df_filtered.empty:
        logging.info(f"Skipping empty DataFrame for {dataset_name}, partition {partition_desc}, peptide {peptide}")
        return  # Skip if there's nothing to process

    logging.info(f"Starting processing for {dataset_name}, partition {partition_desc}, peptide {peptide}. Total rows: {len(df_filtered)}")

    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        data = [(row["peptide_x"], row["tcr_full"]) for _, row in df_subset.iterrows()]
        logging.info(f"Processing batch {start} to {end} for {dataset_name}, partition {partition_desc}, peptide {peptide}")

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            mean_tensor = token_representations[i, 1 : tokens_len - 1].mean(0)
            all_sequence_representations.append(mean_tensor.cpu().numpy())

    output_file_path = os.path.join(
        output_dir,
        f"sequence_representations_{dataset_name}_{partition_desc}_{peptide}.csv",
    )
    
    pd.DataFrame(all_sequence_representations).to_csv(output_file_path, index=False)
    logging.info(f"Saved to {output_file_path}")


# Loop through peptides, partitions, and datasets
for peptide in peptides:
    logging.info(f"\nProcessing peptide: {peptide}")

    for partition in partitions:
        logging.info(f"\nProcessing partition {partition} for binders dataset...")

        # Process binders_only dataset for each partition
        df_binders_only = pd.read_csv(datasets["binders_only"])
        df_filtered_binders_only = df_binders_only[
            (df_binders_only["partition"] == partition)
            & (df_binders_only["peptide_x"] == peptide)
        ]
        
        # Check if filtered DataFrame is empty
        if df_filtered_binders_only.empty:
            logging.info(f"No data found for peptide {peptide}, partition {partition} in binders_only dataset. Skipping...")
        else:
            process_sequences(
                df_filtered_binders_only, peptide, f"partition_{partition}", "binders_only"
            )

        # Process not binders_only dataset for various partition combinations
        logging.info(f"\nProcessing partition combinations for full_sequence dataset...")

        df_full_sequence = pd.read_csv(datasets["full_sequence"])

        # Full partitions (1,2,3,4) except the current one
        partition_combos = [
            (1, 2, 3, 4),
            (0, 2, 3, 4),
            (0, 1, 3, 4),
            (0, 1, 2, 4),
            (0, 1, 2, 3),
        ]
        
        for combo in partition_combos:
            df_filtered_full_sequence = df_full_sequence[
                (df_full_sequence["partition"].isin(combo))
                & (df_full_sequence["peptide_x"] == peptide)
            ]
            partition_desc = f'partitions_{"_".join(map(str, combo))}'

            if df_filtered_full_sequence.empty:
                logging.info(f"No data found for peptide {peptide}, partition combo {combo} in full_sequence dataset. Skipping...")
            else:
                process_sequences(
                    df_filtered_full_sequence, peptide, partition_desc, "full_sequence"
                )

logging.info("\nAll processing complete.")

