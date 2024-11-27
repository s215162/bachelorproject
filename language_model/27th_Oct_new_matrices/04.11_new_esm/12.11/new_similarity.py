import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import sys
import logging

# Set up logging
log_file = "processing_log.txt"
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
)

def calculate_similarity(peptide, weight):
    logging.info(f"Processing peptide: {peptide}, weight: {weight}...")

    # Paths for binders and test files
    train_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{weight.lower()}/sequence_summed_vectors_binders_{peptide}_{weight}_np.csv"
    test_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{weight.lower()}/sequence_summed_vectors_swaps_{peptide}_{weight}_np.csv"

    if not os.path.exists(train_file):
        logging.error(f"Training file not found: {train_file}")
        return
    if not os.path.exists(test_file):
        logging.error(f"Test file not found: {test_file}")
        return

    # Load data
    training_data = pd.read_csv(train_file)
    test_data = pd.read_csv(test_file)

    # Extract embeddings and raw indices
    training_embeddings = training_data.drop(columns=["raw_index"]).values
    test_embeddings = test_data.drop(columns=["raw_index"]).values
    test_raw_indices = test_data["raw_index"].values
    training_raw_indices = training_data["raw_index"].values

    # Normalize embeddings
    training_embeddings_norm = training_embeddings / np.linalg.norm(
        training_embeddings, axis=1, keepdims=True
    )
    test_embeddings_norm = test_embeddings / np.linalg.norm(
        test_embeddings, axis=1, keepdims=True
    )

    # Calculate maximum similarity for each test entry
    results = {}
    for i, test_emb in enumerate(test_embeddings_norm):
        test_raw_index = test_raw_indices[i]

        # Filter training data by excluding identical raw_index
        valid_indices = [
            j for j, train_raw_index in enumerate(training_raw_indices)
            if train_raw_index != test_raw_index
        ]
        filtered_training_embeddings = training_embeddings_norm[valid_indices]

        if len(filtered_training_embeddings) == 0:
            logging.warning(f"No valid comparisons for test sample {test_raw_index}.")
            continue

        similarities = cosine_similarity([test_emb], filtered_training_embeddings)[0]
        max_similarity = similarities.max()
        results[test_raw_index] = max(results.get(test_raw_index, 0), max_similarity)

    # Save results to CSV
    results_df = pd.DataFrame(results.items(), columns=["raw_index", "max_similarity"])
    output_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/cos_sim/{weight.lower()}_{peptide}_similarity_results.csv"
    results_df.to_csv(output_file, index=False)
    logging.info(f"Results saved to {output_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_similarity.py <peptide> <weight>")
        sys.exit(1)

    peptide = sys.argv[1]
    weight = sys.argv[2]
    calculate_similarity(peptide, weight)

