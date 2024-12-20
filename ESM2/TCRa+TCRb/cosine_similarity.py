import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import logging
import sys
import os

def setup_logging(peptide):
    """Configure logging to track similarity calculations for specific peptides"""
    # change path baseed on user
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/new_similarity_processing_{peptide}.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    # Add console output for monitoring
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def calculate_similarities(test_file, train_file, output_file):
    """Calculate cosine similarities between ESM embeddings.
    
    Compares each test sequence against binder-only training sequences
    within the same partition. Uses summed embeddings for comparison.
    """
    # Load test and training data from NPZ files
    test_data = np.load(test_file, allow_pickle=True)
    train_data = np.load(train_file, allow_pickle=True)
    
    # Extract arrays and metadata
    test_embeddings = test_data['embeddings']
    test_raw_indexes = test_data['raw_indexes']
    test_binder_values = test_data['binder_values']
    test_partition_values = test_data['partition_values']
    
    train_embeddings = train_data['embeddings']
    train_raw_indexes = train_data['raw_indexes']
    train_partition_values = train_data['partition_values']
    
    results = []
    
    # Convert embeddings matrices to vectors
    test_summed = [np.sum(emb, axis=0) for emb in test_embeddings]
    train_summed = [np.sum(emb, axis=0) for emb in train_embeddings]
    
    # Compare each test sequence against all training sequences
    for i, test_vec in enumerate(test_summed):
        test_raw_index = test_raw_indexes[i]
        test_binder = test_binder_values[i]
        test_partition = test_partition_values[i]
        
        max_similarity = -1
        
        # Find most similar training sequence
        for j, train_vec in enumerate(train_summed):
            train_raw_index = train_raw_indexes[j]
            train_partition = train_partition_values[j]
            
            # Only compare within same partition and different sequences
            if train_partition != test_partition or train_raw_index == test_raw_index:
                continue
                
            # Calculate cosine similarity between vectors
            similarity = cosine_similarity(
                test_vec.reshape(1, -1),
                train_vec.reshape(1, -1)
            )[0][0]
            
            if similarity > max_similarity:
                max_similarity = similarity
        
        # Store results for this test sequence
        results.append({
            'raw_index': test_raw_index,
            'max_similarity': max_similarity,
            'binder': test_binder,
            'partition': test_partition
        })
    
    # Save similarity results
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False)
    logging.info(f"Saved results to {output_file}")

def process_cdr3_similarities(peptide):
    """Process similarities for CDR3 regions independently.
    
    Analyzes alpha and beta chain CDR3 regions separately to assess
    their individual contributions to binding prediction.
    """
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains/confirming"
    os.makedirs(base_dir, exist_ok=True)
    
    # Process each CDR3 region independently
    for region in ['CDR3a', 'CDR3b']:
        logging.info(f"Processing {region}")
        
        test_file = os.path.join(base_dir, f"{peptide}_{region}_full.npz")
        train_file = os.path.join(base_dir, f"{peptide}_{region}_binder.npz")
        output_file = os.path.join(base_dir, f"{peptide}_{region}_similarities.csv")
        
        if not (os.path.exists(test_file) and os.path.exists(train_file)):
            logging.error(f"Missing required files for {region}")
            continue
            
        try:
            calculate_similarities(test_file, train_file, output_file)
            logging.info(f"Completed processing for {region}")
        except Exception as e:
            logging.error(f"Error processing {region}: {e}")

def calculate_combined_cdr3_similarities(peptide):
    """Calculate combined CDR3 region similarities.
    
    Assesses combined predictive power of alpha and beta CDR3 regions
    using both sum and mean-based approaches for vector combination.
    """
    # change path based on user
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains/confirming"
    
    # Load CDR3 data for both chains
    cdr3a_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_full.npz"), allow_pickle=True)
    cdr3b_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_full.npz"), allow_pickle=True)
    cdr3a_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_binder.npz"), allow_pickle=True)
    cdr3b_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_binder.npz"), allow_pickle=True)
    
    # Get partition information for comparison control
    test_partition_values = cdr3a_test['partition_values']
    train_partition_values = cdr3a_train['partition_values']
    
    # Process embeddings using both sum and mean approaches
    cdr3a_test_summed = [np.sum(emb, axis=0) for emb in cdr3a_test['embeddings']]
    cdr3b_test_summed = [np.sum(emb, axis=0) for emb in cdr3b_test['embeddings']]
    cdr3a_train_summed = [np.sum(emb, axis=0) for emb in cdr3a_train['embeddings']]
    cdr3b_train_summed = [np.sum(emb, axis=0) for emb in cdr3b_train['embeddings']]
    
    cdr3a_test_mean = [np.mean(emb, axis=0) for emb in cdr3a_test['embeddings']]
    cdr3b_test_mean = [np.mean(emb, axis=0) for emb in cdr3b_test['embeddings']]
    cdr3a_train_mean = [np.mean(emb, axis=0) for emb in cdr3a_train['embeddings']]
    cdr3b_train_mean = [np.mean(emb, axis=0) for emb in cdr3b_train['embeddings']]
    
    results_sum = []
    results_mean = []
    
    # Process each test sequence
    for i in range(len(cdr3a_test_summed)):
        test_raw_index = cdr3a_test['raw_indexes'][i]
        test_binder = cdr3a_test['binder_values'][i]
        test_partition = test_partition_values[i]
        
        # Track best similarities for both approaches
        best_combined_score_sum = -2
        best_similarity_a_sum = -1
        best_similarity_b_sum = -1
        
        best_combined_score_mean = -2
        best_similarity_a_mean = -1
        best_similarity_b_mean = -1
        
        # Compare against all training sequences
        for j in range(len(cdr3a_train_summed)):
            train_raw_index = cdr3a_train['raw_indexes'][j]
            train_partition = train_partition_values[j]
            
            if train_partition != test_partition or train_raw_index == test_raw_index:
                continue
            
            # calculate similarities by summing over the columns
            similarity_a_sum = cosine_similarity(
                cdr3a_test_summed[i].reshape(1, -1),
                cdr3a_train_summed[j].reshape(1, -1)
            )[0][0]
            
            similarity_b_sum = cosine_similarity(
                cdr3b_test_summed[i].reshape(1, -1),
                cdr3b_train_summed[j].reshape(1, -1)
            )[0][0]
            
            combined_score_sum = similarity_a_sum + similarity_b_sum
            
            # calculate similarities by averaging the columns
            similarity_a_mean = cosine_similarity(
                cdr3a_test_mean[i].reshape(1, -1),
                cdr3a_train_mean[j].reshape(1, -1)
            )[0][0]
            
            similarity_b_mean = cosine_similarity(
                cdr3b_test_mean[i].reshape(1, -1),
                cdr3b_train_mean[j].reshape(1, -1)
            )[0][0]
            
            combined_score_mean = similarity_a_mean + similarity_b_mean
            
            # Update best scores if needed
            if combined_score_sum > best_combined_score_sum:
                best_combined_score_sum = combined_score_sum
                best_similarity_a_sum = similarity_a_sum
                best_similarity_b_sum = similarity_b_sum
                
            if combined_score_mean > best_combined_score_mean:
                best_combined_score_mean = combined_score_mean
                best_similarity_a_mean = similarity_a_mean
                best_similarity_b_mean = similarity_b_mean
        
        # Store resultss
        results_sum.append({
            'raw_index': test_raw_index,
            'max_similarity_a': best_similarity_a_sum,
            'max_similarity_b': best_similarity_b_sum,
            'sum_max_similarity': best_combined_score_sum,
            'binder': test_binder,
            'partition': test_partition
        })
        
        results_mean.append({
            'raw_index': test_raw_index,
            'max_similarity_a': best_similarity_a_mean,
            'max_similarity_b': best_similarity_b_mean,
            'sum_max_similarity': best_combined_score_mean,
            'binder': test_binder,
            'partition': test_partition
        })

    # Save both sets 
    pd.DataFrame(results_sum).to_csv(
        os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_sum.csv"), 
        index=False
    )
    pd.DataFrame(results_mean).to_csv(
        os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_mean.csv"), 
        index=False
    )
    logging.info("Saved both sum and mean CDR3 results")

def calculate_all_cdr_similarities(peptide):
    """Calculate similarities using all CDR regions.
    
    Combines all CDR regions (1-3, alpha and beta) using both weighted
    and unweighted approaches. CDR3 regions receive 4x weight in the
    weighted approach.
    """
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains/confirming"
    cdr_regions = ['CDR1a', 'CDR1b', 'CDR2a', 'CDR2b', 'CDR3a', 'CDR3b']
    
    # Load embeddings
    test_data = {}
    train_data = {}
    for region in cdr_regions:
        test_path = os.path.join(base_dir, f"{peptide}_{region}_full.npz")
        train_path = os.path.join(base_dir, f"{peptide}_{region}_binder.npz")
        
        if not (os.path.exists(test_path) and os.path.exists(train_path)):
            logging.error(f"Missing required files for {region}")
            return
            
        test_data[region] = np.load(test_path, allow_pickle=True)
        train_data[region] = np.load(train_path, allow_pickle=True)
    
    # Process embeddings
    test_summed = {}
    train_summed = {}
    for region in cdr_regions:
        test_summed[region] = [np.sum(emb, axis=0) for emb in test_data[region]['embeddings']]
        train_summed[region] = [np.sum(emb, axis=0) for emb in train_data[region]['embeddings']]
    
    results = []
    
    # Process each test sequence
    for i in range(len(test_summed['CDR1a'])):
        test_raw_index = test_data['CDR1a']['raw_indexes'][i]
        test_binder = test_data['CDR1a']['binder_values'][i]
        test_partition = test_data['CDR1a']['partition_values'][i]
        
        best_weighted_score = float('-inf')
        best_similarities = {region: -1 for region in cdr_regions}
        best_unweighted_sum = -1
        
        # Compare
        for j in range(len(train_summed['CDR1a'])):
            train_raw_index = train_data['CDR1a']['raw_indexes'][j]
            train_partition = train_data['CDR1a']['partition_values'][j]
            
            if train_partition != test_partition or train_raw_index == test_raw_index:
                continue
            
            # Calculate similarities for regions
            current_similarities = {}
            for region in cdr_regions:
                similarity = cosine_similarity(
                    test_summed[region][i].reshape(1, -1),
                    train_summed[region][j].reshape(1, -1)
                )[0][0]
                current_similarities[region] = similarity
            
            # Calculate combined scores
            unweighted_sum = sum(current_similarities.values())
            
            # 4x weight for CDR3 regions
            weighted_sum = (
                current_similarities['CDR1a'] + 
                current_similarities['CDR1b'] + 
                current_similarities['CDR2a'] + 
                current_similarities['CDR2b'] + 
                4 * current_similarities['CDR3a'] + 
                4 * current_similarities['CDR3b']
            )
            
            # Update scores
            if weighted_sum > best_weighted_score:
                best_weighted_score = weighted_sum
                best_similarities = current_similarities.copy()
                best_unweighted_sum = unweighted_sum
        
        # Store results
        result = {
            'raw_index': test_raw_index,
            'binder': test_binder,
            'weighted_sum': best_weighted_score,
            'unweighted_sum': best_unweighted_sum,
            'partition': test_partition
        }
        # Add individual region similarities to results
        for region in cdr_regions:
            result[f'max_similarity_{region}'] = best_similarities[region]
        
        results.append(result)
    
    # Save combined results
    output_file = os.path.join(base_dir, f"{peptide}_all_CDR_combined_similarities.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved all CDR combined results to {output_file}")

def main():
    """Process TCR similarity analysis for a specific peptide.
    
    Performs three levels of analysis:
    1. Individual CDR3 regions
    2. Combined CDR3 (sum and mean based)
    3. All CDR regions together (weighted/unweighted)
    """
    # Verify command line arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptide = sys.argv[1]
    setup_logging(peptide)
    
    try:
        # Run fulll analysis
        process_cdr3_similarities(peptide)
        logging.info("Completed individual CDR3 analysis")
        
        calculate_combined_cdr3_similarities(peptide)
        logging.info("Completed combined CDR3 analysis")
        
        calculate_all_cdr_similarities(peptide)
        logging.info("Completed all CDR combined analysis")
        
        logging.info(f"All processing completed for peptide {peptide}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()