import torch
import esm
import pandas as pd
import os
import logging

# Setup logging
log_file = ("/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/logs/alpha_full_esm_processing.log")

logging.basicConfig(
    filename=log_file,
    filemode="a",  # Append mode
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)
console = logging.StreamHandler()  # Log to console as well
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

# Define peptides and partitions
peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
partitions = [0, 1, 2, 3, 4]

# Define the datasets
datasets = {
    'binders': '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv',
    'swaps': '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv'
}

# Define output directory
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/full"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Load datasets only once
df_binders = pd.read_csv(datasets['binders'])
swaps = pd.read_csv(datasets['swaps'])

# Function to process and save sequence representations
def process_sequences(df_filtered, peptide, partition_desc, dataset_name):
    all_sequence_representations = []
    raw_indexes = []
    batch_size = 25  # Adjust based on memory capacity
    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        data = [(row["peptide_x"], row["TCRa"]) for _, row in df_subset.iterrows()]
        raw_indexes.extend(df_subset["raw_index"].tolist())  # Collect raw_index values
        logging.info(f"Processing rows {start} to {end} for {dataset_name}, partition {partition_desc}, peptide {peptide}")

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            mean_tensor = token_representations[i, 1:tokens_len - 1].mean(0)
            all_sequence_representations.append(mean_tensor.cpu().numpy())

    # Create a DataFrame with raw_index and sequence representations
    output_df = pd.DataFrame(all_sequence_representations)
    output_df.insert(0, "raw_index", raw_indexes)  # Insert raw_index as the first column

    output_file_path = os.path.join(output_dir, f'sequence_representations_{dataset_name}_{partition_desc}_{peptide}.csv')
    output_df.to_csv(output_file_path, index=False)
    logging.info(f'Saved to {output_file_path}')

# Loop through peptides and partitions
for peptide in peptides:
    for partition in partitions:
        # Process binders dataset for each partition
        df_filtered_binders = df_binders[(df_binders['partition'] == partition) & (df_binders['peptide_x'] == peptide)]
        logging.info(f"Processing binders for partition {partition} and peptide {peptide}")
        process_sequences(df_filtered_binders, peptide, f'partition_{partition}', 'binders')

        # Process swaps dataset for various partition combinations
        # Full partitions (1,2,3,4) except the current one
    partition_combos = [(1, 2, 3, 4), (0, 2, 3, 4), (0, 1, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)]
    for combo in partition_combos:
            df_filtered_swaps = swaps[(swaps['partition'].isin(combo)) & (swaps['peptide_x'] == peptide)]
            partition_desc = f'partitions_{"_".join(map(str, combo))}'
            logging.info(f"Processing swaps for partitions {partition_desc} and peptide {peptide}")
            process_sequences(df_filtered_swaps, peptide, partition_desc, 'swap')

