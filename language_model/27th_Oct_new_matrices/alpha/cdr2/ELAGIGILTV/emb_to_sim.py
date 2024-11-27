import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import matplotlib.pyplot as plt

peptide = "ELAGIGILTV"
type = "alpha"
sequence = "cdr2"
output_dis = "cdr2a"
partitions = [
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_0_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_1_2_3_4_"
        + peptide
        + ".csv",
    },
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_1_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_0_2_3_4_"
        + peptide
        + ".csv",
    },
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_2_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_0_1_3_4_"
        + peptide
        + ".csv",
    },
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_3_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_0_1_2_4_"
        + peptide
        + ".csv",
    },
    {
        "train": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_binders_partition_4_"
        + peptide
        + ".csv",
        "test": "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/"
        + type
        + "/"
        + sequence
        + "/"
        + peptide
        + "/sequence_representations_swap_partitions_0_1_2_3_"
        + peptide
        + ".csv",
    },
]

# Create directory to save histograms if it doesn't exist
output_dir = "similarity_histograms"
os.makedirs(output_dir, exist_ok=True)

# Initialize dictionary to store max similarity results
results = {}

for round_num, round in enumerate(partitions):

    # Load the training and test embeddings
    training_data = pd.read_csv(round["train"])
    test_data = pd.read_csv(round["test"])

    # Extract and normalize embeddings
    training_embeddings = training_data.drop(columns=["raw_index"]).values
    test_embeddings = test_data.drop(columns=["raw_index"]).values
    test_raw_indices = test_data["raw_index"].values

    # Normalize embeddings for cosine similarity
    training_embeddings_norm = training_embeddings / np.linalg.norm(
        training_embeddings, axis=1, keepdims=True
    )
    test_embeddings_norm = test_embeddings / np.linalg.norm(
        test_embeddings, axis=1, keepdims=True
    )

    # Calculate similarities and find the maximum similarity for each test point
    all_similarities = []
    for i, test_emb in enumerate(test_embeddings_norm):
        similarities = cosine_similarity([test_emb], training_embeddings_norm)[0]
        max_similarity = similarities.max()
        raw_index = test_raw_indices[i]

        # Store max similarity for each raw_index
        if raw_index in results:
            results[raw_index] = max(results[raw_index], max_similarity)
        else:
            results[raw_index] = max_similarity

        # Add similarities to the list for distribution analysis
        all_similarities.extend(similarities)

    # Plot the distribution of similarities for this round and save as PNG
    plt.figure(figsize=(8, 6))
    plt.hist(all_similarities, bins=50, alpha=0.7, color="blue", edgecolor="black")
    plt.title(f"Cosine Similarity Distribution - Round {round_num + 1}")
    plt.xlabel("Cosine Similarity")
    plt.ylabel("Frequency")
    plt.grid(True)

    # Save histogram
    plt.savefig(
        os.path.join(output_dir, f"similarity_distribution_round_{round_num + 1}.png")
    )
    plt.close()

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(
    list(results.items()), columns=["raw_index", "max_similarity"]
)
results_df.to_csv(
    "max_cosine_similarity_results_" + peptide + "_" + output_dis + ".csv", index=False
)

