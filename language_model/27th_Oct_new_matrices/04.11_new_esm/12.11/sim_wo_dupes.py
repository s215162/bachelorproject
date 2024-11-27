import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import logging
from datetime import datetime
from multiprocessing import Pool, cpu_count

# Set up logging
log_file = '/net/mimer/mnt/tank/projects2/emison/language_model/processing_log.txt'
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
)

# Peptides and weights
peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
weights = ["CDR1a", "CDR1b", "CDR2a", "CDR2b", "CDR3a", "CDR3b", "TCRa", "TCRb"]

# Load the full sequence data once
logging.info("Loading full sequence data...")
full_sequence_df = pd.read_csv('/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv')
logging.info("Full sequence data loaded successfully.")

# Processing function for each peptide-weight pair
def process_peptide_weight(args):
    peptide, weight = args
    logging.info(f"Processing peptide: {peptide}, weight: {weight}...")

    # File paths
    train_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/{weight.lower()}/sequence_summed_vectors_binders_{peptide}_{weight}_np.csv"
    test_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/{weight.lower()}/sequence_summed_vectors_swaps_{peptide}_{weight}_np.csv"

    # Check if files exist
    if not os.path.exists(train_file) or not os.path.exists(test_file):
        logging.error(f"Missing file for peptide: {peptide}, weight: {weight}")
        return

    # Load data
    training_data = pd.read_csv(train_file)
    test_data = pd.read_csv(test_file)

    training_embeddings = training_data.drop(columns=["raw_index"]).values
    test_embeddings = test_data.drop(columns=["raw_index"]).values
    test_raw_indices = test_data["raw_index"].values
    training_raw_indices = training_data["raw_index"].values

    # Normalize embeddings
    training_embeddings /= np.linalg.norm(training_embeddings, axis=1, keepdims=True)
    test_embeddings /= np.linalg.norm(test_embeddings, axis=1, keepdims=True)

    # Compute max similarities
    results = {}
    for i, test_emb in enumerate(test_embeddings):
        test_raw_index = test_raw_indices[i]

        # Filter by raw_index
        valid_indices = [j for j, train_raw_index in enumerate(training_raw_indices)
                         if train_raw_index != test_raw_index]
        if not valid_indices:
            continue

        filtered_training_embeddings = training_embeddings[valid_indices]
        similarities = cosine_similarity([test_emb], filtered_training_embeddings)[0]
        max_similarity = similarities.max()
        results[test_raw_index] = max(results.get(test_raw_index, 0), max_similarity)

    # Convert to DataFrame
    results_df = pd.DataFrame(list(results.items()), columns=["raw_index", "max_similarity"])

    # Merge with full_sequence_df to include binders
    merged_df = pd.merge(
        results_df,
        full_sequence_df[['raw_index', 'peptide_x', 'binders']],
        how='left',
        on='raw_index'
    ).drop(columns=['peptide_x'])

    # Save results
    output_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/{weight.lower()}/max_cosine_similarity_results_{peptide}_{weight}.csv"
    merged_df.to_csv(output_file, index=False)
    logging.info(f"Results for peptide: {peptide}, weight: {weight} saved to {output_file}.")

# Main function
if __name__ == "__main__":
    tasks = [(peptide, weight) for peptide in peptides for weight in weights]
    with Pool(cpu_count()) as pool:
        pool.map(process_peptide_weight, tasks)
    logging.info("All processing complete.")

