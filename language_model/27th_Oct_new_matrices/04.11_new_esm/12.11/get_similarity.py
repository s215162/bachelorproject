# get_similarity

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import os
import logging
from datetime import datetime

# Set up logging
log_file = '/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/processing_log.txt'
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
)

# List of peptides and weights
peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
weights = ["CDR1a", "CDR1b", "CDR2a", "CDR2b", "CDR3a", "CDR3b", "TCRa", "TCRb"]

# Load the full sequence data once (to avoid loading it multiple times)
logging.info("Loading full sequence data...")
full_sequence_df = pd.read_csv('/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv')
logging.info("Full sequence data loaded successfully.")

# Loop through peptides and weights
for peptide in peptides:
    for weight in weights:
        logging.info(f"Processing peptide: {peptide}, weight: {weight}...")

        # Construct paths for binders and data files
        train_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{weight.lower()}/sequence_summed_vectors_binders_{peptide}_{weight}_np.csv"
        test_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{weight.lower()}/sequence_summed_vectors_swaps_{peptide}_{weight}_np.csv"

        # Check if training and test files exist
        if not os.path.exists(train_file):
            logging.error(f"Training file not found: {train_file}")
            continue
        if not os.path.exists(test_file):
            logging.error(f"Test file not found: {test_file}")
            continue

        logging.info(f"Loading training data from {train_file}...")
        training_data = pd.read_csv(train_file)
        logging.info("Training data loaded successfully.")

        logging.info(f"Loading test data from {test_file}...")
        test_data = pd.read_csv(test_file)
        logging.info("Test data loaded successfully.")

        # Extract embeddings and raw indices
        training_embeddings = training_data.drop(columns=["raw_index"]).values
        test_embeddings = test_data.drop(columns=["raw_index"]).values
        test_raw_indices = test_data["raw_index"].values

        # Normalize embeddings
        logging.info("Normalizing embeddings...")
        training_embeddings_norm = training_embeddings / np.linalg.norm(
            training_embeddings, axis=1, keepdims=True
        )
        test_embeddings_norm = test_embeddings / np.linalg.norm(
            test_embeddings, axis=1, keepdims=True
        )
              
        # Calculate similarities and save max similarity for each test entry
        results = {}
        for i, test_emb in enumerate(test_embeddings_norm):
            similarities = cosine_similarity([test_emb], training_embeddings_norm)[0]
            max_similarity = similarities.max()
            raw_index = test_raw_indices[i]

            # Store maximum similarity for each raw index
            results[raw_index] = max(results.get(raw_index, 0), max_similarity)

        # Convert results to DataFrame
        results_df = pd.DataFrame(
            list(results.items()), columns=["raw_index", "max_similarity"]
        )
        logging.info(f"Max cosine similarities calculated for {len(results)} entries.")

        # Merge with full_sequence_df to get the 'binders' column
        logging.info(f"Merging results with full sequence data for peptide: {peptide}...")
        merged_df = pd.merge(
            results_df,
            full_sequence_df[['raw_index', 'peptide_x', 'binders']],
            how='left',
            left_on=['raw_index', 'peptide'],
            right_on=['raw_index', 'peptide_x']
        )

        # Drop 'peptide_x' column after merge as it is not needed
        merged_df = merged_df.drop(columns=['peptide_x'])

        # Save the final result to CSV
        output_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{weight.lower()}/max_cosine_similarity_results_{peptide}_{weight}.csv"
        merged_df.to_csv(output_file, index=False)
        logging.info(f"Processed {peptide} with weight {weight}. Results saved to {output_file}.")

logging.info("Processing complete.")
