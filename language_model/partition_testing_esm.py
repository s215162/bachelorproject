import torch
import esm
import pandas as pd

# Load the dataset
datafile_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data_binders_only.csv"  # Path to your data file
df = pd.read_csv(datafile_path)

# Filter for rows where partition == 0 and peptide_x == "ELAGIGILTV"
df_filtered = df[(df['partition'] != 1) & (df['peptide_x'] == "ELAGIGILTV")]

# Check how many rows will be processed
print(f"Number of rows to process: {len(df_filtered)}")

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Initialize a list to store all sequence representations
all_sequence_representations = []

# Define the batch size
batch_size = 25  # Adjust based on memory capacity

# Process the filtered DataFrame in smaller batches
for start in range(0, len(df_filtered), batch_size):
    end = min(start + batch_size, len(df_filtered))
    df_subset = df_filtered.iloc[start:end]
    
    # Prepare the data from your file: tuples of (peptide, tcr_full)
    data = [(row["peptide_x"], row["tcr_full"]) for _, row in df_subset.iterrows()]
    print(f"Processing rows {start} to {end}")

    # Convert the data to the format needed by the model
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)

    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    for i, tokens_len in enumerate(batch_lens):
        # Convert tensor to NumPy array
        mean_tensor = token_representations[i, 1:tokens_len - 1].mean(0)
        all_sequence_representations.append(mean_tensor.cpu().numpy())  # Move to CPU and convert to NumPy

# Save all representations to a new CSV file
output_file_path = "/net/mimer/mnt/tank/projects2/emison/language_model/sequence_representations_BINDERS_partition_!1_ELAGIGILTV.csv"
pd.DataFrame(all_sequence_representations).to_csv(output_file_path, index=False)

print(f'Sequence representations for partitions 0, 2,3, 4 and peptide_x "ELAGIGILTV" saved to {output_file_path}')

