import torch
import esm
import pandas as pd

# Load the dataset
datafile_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"  # Path to your data file
df = pd.read_csv(datafile_path)

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Prepare the data from your file: tuples of (peptide, tcr_full)
# You can use either 'peptide_x' or 'peptide_y' as the name, depending on your use case.
data = [(row["peptide_x"], row["tcr_full"]) for _, row in df.iterrows()]
print("Data retrieved")

# Convert the data to the format needed by the model
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)
token_representations = results["representations"][33]

# Generate per-sequence representations via averaging
# NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
sequence_representations = []
for i, tokens_len in enumerate(batch_lens):
    sequence_representations.append(
        token_representations[i, 1:tokens_len - 1].mean(0).numpy()  # Convert to numpy array for saving
    )

# Prepare the output DataFrame
output_df = pd.DataFrame(sequence_representations)
output_df.columns = [f"representation_{i+1}" for i in range(output_df.shape[1])]  # Name the columns
output_df['peptide_x'] = [item[0] for item in data]  # Add peptide_x for reference
output_df['tcr_full'] = [item[1] for item in data]  # Add tcr_full for reference

# Save the output DataFrame to a CSV file
output_file_path = "/net/mimer/mnt/tank/projects2/emison/language_model/sequence_representations.csv"
output_df.to_csv(output_file_path, index=False)

print(f"Sequence representations saved to {output_file_path}")

