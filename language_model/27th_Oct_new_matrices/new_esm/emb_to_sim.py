import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Generate similarity matrix based on embeddings.")
parser.add_argument("--peptide", type=str, required=True, help="Specify the peptide sequence.")
parser.add_argument("--type", type=str, required=True, help="Specify the type (e.g., alpha or beta).")
parser.add_argument("--sequence", type=str, required=True, help="Specify the sequence part (e.g., cdr1, cdr2, etc.).")
parser.add_argument("--output_dis", type=str, required=True, help="Specify the output file prefix or name.")
parser.add_argument("--output_dir", type=str, default="output", help="Specify the base output directory.")

# Parse arguments
args = parser.parse_args()

# Assign variables from parsed arguments
peptide = args.peptide
type_ = args.type  # Change 'type' to 'type_'
sequence = args.sequence
output_dis = args.output_dis
base_output_dir = args.output_dir

# Construct output directory path based on peptide
output_dir = os.path.join(base_output_dir, peptide)
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Define partitions with corrected 'type_' usage
partitions = [
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type_
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_0_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type_
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_1_2_3_4_"
        + peptide
        + ".csv",
    },
    # Continue for other partitions...
]

# Create directory to save histograms if it doesn't exist
output_dir = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
    + type_
    + "/"
    + sequence
    + "/"
    + peptide
    + "/similarity_histograms_"
    + peptide
    + "_"
    + output_dis
)
os.makedirs(output_dir, exist_ok=True)

# Remaining code
results = {}
for round_num, round in enumerate(partitions):
    training_data = pd.read_csv(round["train"])
    test_data = pd.read_csv(round["test"])
    training_embeddings = training_data.drop(columns=["raw_index"]).values
    test_embeddings = test_data.drop(columns=["raw_index"]).values
    test_raw_indices = test_data["raw_index"].values

    training_embeddings_norm = training_embeddings / np.linalg.norm(training_embeddings, axis=1, keepdims=True)
    test_embeddings_norm = test_embeddings / np.linalg.norm(test_embeddings, axis=1, keepdims=True)

    all_similarities = []
    for i, test_emb in enumerate(test_embeddings_norm):
        similarities = cosine_similarity([test_emb], training_embeddings_norm)[0]
        max_similarity = similarities.max()
        raw_index = test_raw_indices[i]

        if raw_index in results:
            results[raw_index] = max(results[raw_index], max_similarity)
        else:
            results[raw_index] = max_similarity

        all_similarities.extend(similarities)

    plt.figure(figsize=(8, 6))
    plt.hist(all_similarities, bins=50, alpha=0.7, color="blue", edgecolor="black")
    plt.title(f"Cosine Similarity Distribution - Round {round_num + 1}")
    plt.xlabel("Cosine Similarity")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"similarity_distribution_round_{round_num + 1}.png"))
    plt.close()

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(list(results.items()), columns=["raw_index", "max_similarity"])
results_df.to_csv(os.path.join(base_output_dir, f"max_cosine_similarity_results_{peptide}_{output_dis}.csv"), index=False)

