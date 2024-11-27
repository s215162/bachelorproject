import torch
import esm
import pandas as pd

# Load the dataset
datafile_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"  # Path to your data file
df = pd.read_csv(datafile_path)

# Select the first 10 rows
df_subset = df.head(10)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()  # Change to the desired model
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Prepare the data from your file: tuples of (peptide, tcr_full)
data = [(row["peptide_x"], row["tcr_full"]) for _, row in df_subset.iterrows()]
print("Data retrieved")

# Convert the data to the format needed by the model
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)

token_representations = results["representations"][33]

# Generate per-sequence representations via averaging
sequence_representations = []
for i, tokens_len in enumerate(batch_lens):
    # Convert tensor to NumPy array
    mean_tensor = token_representations[i, 1:tokens_len - 1].mean(0)
    sequence_representations.append(mean_tensor.cpu().numpy())  # Move to CPU and convert to NumPy

# Save the representations to a new CSV file
output_file_path = "/net/mimer/mnt/tank/projects2/emison/language_model/sequence_representations.csv"
pd.DataFrame(sequence_representations).to_csv(output_file_path, index=False)

print(f'Sequence representations saved to {output_file_path}')

