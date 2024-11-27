import torch
import esm
import pandas as pd
import os

# Define peptides and partitions
peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
partitions = [0, 1, 2, 3, 4]

# Define the datasets
datasets = {
    'binders_only': '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv',
    'full_sequence': '/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv'
}

# Define output directory
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/esm_unweighed"
os.makedirs(output_dir, exist_ok=True)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Function to process and save sequence representations
def process_sequences(df_filtered, peptide, partition_desc, dataset_name):
    all_sequence_representations = []
    batch_size = 25  # Adjust based on memory capacity

    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        df_subset = df_filtered.iloc[start:end]

        data = [(row["peptide_x"], row["tcr_full"]) for _, row in df_subset.iterrows()]
        print(f"Processing rows {start} to {end} for {dataset_name}, partition {partition_desc}, peptide {peptide}")

        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)

        token_representations = results["representations"][33]

        for i, tokens_len in enumerate(batch_lens):
            mean_tensor = token_representations[i, 1:tokens_len - 1].mean(0)
            all_sequence_representations.append(mean_tensor.cpu().numpy())

    output_file_path = os.path.join(output_dir, f'sequence_representations_{dataset_name}_{partition_desc}_{peptide}.csv')
    pd.DataFrame(all_sequence_representations).to_csv(output_file_path, index=False)
    print(f'Saved to {output_file_path}')

# Loop through peptides, partitions, and datasets
for peptide in peptides:
    for partition in partitions:
        # Process binders_only dataset for each partition
        df_binders_only = pd.read_csv(datasets['binders_only'])
        df_filtered_binders_only = df_binders_only[(df_binders_only['partition'] == partition) & (df_binders_only['peptide_x'] == peptide)]
        process_sequences(df_filtered_binders_only, peptide, f'partition_{partition}', 'binders_only')

        # Process not binders_only dataset for various partition combinations
        df_full_sequence = pd.read_csv(datasets['full_sequence'])

        # Full partitions (1,2,3,4) except the current one
        partition_combos = [(1, 2, 3, 4), (0, 2, 3, 4), (0, 1, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)]
        for combo in partition_combos:
            df_filtered_full_sequence = df_full_sequence[(df_full_sequence['partition'].isin(combo)) & (df_full_sequence['peptide_x'] == peptide)]
            partition_desc = f'partitions_{"_".join(map(str, combo))}'
            process_sequences(df_filtered_full_sequence, peptide, partition_desc, 'full_sequence')

