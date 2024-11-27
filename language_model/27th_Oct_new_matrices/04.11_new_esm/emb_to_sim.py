import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Generate similarity matrix based on embeddings."
)
parser.add_argument(
    "--peptide", type=str, required=True, help="Specify the peptide sequence."
)
parser.add_argument(
    "--type", type=str, required=True, help="Specify the type (e.g., alpha or beta)."
)
parser.add_argument(
    "--sequence",
    type=str,
    required=True,
    help="Specify the sequence part (e.g., cdr1, cdr2, etc.).",
)
parser.add_argument(
    "--output_dis",
    type=str,
    required=True,
    help="Specify the output file prefix or name.",
)
parser.add_argument(
    "--output_dir",
    type=str,
    default="output",
    help="Specify the base output directory.",
)

# Parse arguments
args = parser.parse_args()

# Assign variables from parsed arguments
peptide = args.peptide
type_ = args.type  # Change 'type' to 'type_'
sequence = args.sequence
output_dis = args.output_dis
base_output_dir = args.output_dir

# Construct paths for binders and data files
train_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/{sequence}/sequence_embeddings_binders_{peptide}.csv"
test_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/{sequence}/sequence_embeddings_data_{peptide}.csv"

# Create output directory for similarity histograms
output_hist_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/{type_}/{sequence}/similarity_histograms_{peptide}_{output_dis}"
os.makedirs(output_hist_dir, exist_ok=True)

# Load data
training_data = pd.read_csv(train_file)
test_data = pd.read_csv(test_file)

# Extract embeddings and raw indices
training_embeddings = training_data.drop(columns=["raw_index"]).values
test_embeddings = test_data.drop(columns=["raw_index"]).values
test_raw_indices = test_data["raw_index"].values

# Normalize embeddings
training_embeddings_norm = training_embeddings / np.linalg.norm(
    training_embeddings, axis=1, keepdims=True
)
test_embeddings_norm = test_embeddings / np.linalg.norm(
    test_embeddings, axis=1, keepdims=True
)

# Calculate similarities and save max similarity for each test entry
results = {}
all_similarities = []
for i, test_emb in enumerate(test_embeddings_norm):
    similarities = cosine_similarity([test_emb], training_embeddings_norm)[0]
    max_similarity = similarities.max()
    raw_index = test_raw_indices[i]

    # Store maximum similarity for each raw index
    results[raw_index] = max(results.get(raw_index, 0), max_similarity)
    all_similarities.extend(similarities)

# Plot and save histogram of cosine similarities
plt.figure(figsize=(8, 6))
plt.hist(all_similarities, bins=50, alpha=0.7, color="blue", edgecolor="black")
plt.title("Cosine Similarity Distribution")
plt.xlabel("Cosine Similarity")
plt.ylabel("Frequency")
plt.grid(True)
plt.savefig(os.path.join(output_hist_dir, f"similarity_distribution.png"))
plt.close()

# Save results to CSV
results_df = pd.DataFrame(
    list(results.items()), columns=["raw_index", "max_similarity"]
)
results_df.to_csv(
    os.path.join(
        base_output_dir, f"max_cosine_similarity_results_{peptide}_{output_dis}.csv"
    ),
    index=False,
)

