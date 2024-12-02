import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import logging
import sys
import os

def setup_logging(peptide):
    """Set up logging configuration."""
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/similarity_processing_{peptide}.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

# TASK 1: Individual CDR3 Analysis
################################

def calculate_similarities(test_file, train_file, output_file):
    """Calculate cosine similarities between test and train arrays."""
    # Load data
    test_data = np.load(test_file, allow_pickle=True)
    train_data = np.load(train_file, allow_pickle=True)
    
    test_embeddings = test_data['embeddings']
    test_raw_indexes = test_data['raw_indexes']
    test_binder_values = test_data['binder_values']
    
    train_embeddings = train_data['embeddings']
    train_raw_indexes = train_data['raw_indexes']
    
    results = []
    
    for i, test_embedding in enumerate(test_embeddings):
        test_raw_index = test_raw_indexes[i]
        test_binder = test_binder_values[i]
        
        max_similarity = -1
        
        for j, train_embedding in enumerate(train_embeddings):
            train_raw_index = train_raw_indexes[j]
            
            if train_raw_index == test_raw_index:
                continue
                
            similarity = cosine_similarity(
                test_embedding.reshape(1, -1),
                train_embedding.reshape(1, -1)
            )[0][0]
            
            if similarity > max_similarity:
                max_similarity = similarity
        
        results.append({
            'raw_index': test_raw_index,
            'max_similarity': max_similarity,
            'binder': test_binder
        })
    
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False)
    logging.info(f"Saved results to {output_file}")

def process_cdr3_similarities(peptide):
    """Process similarities for both CDR3a and CDR3b individually."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    
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

# TASK 2: Combined CDR3 Analysis
##############################

def calculate_combined_cdr3_similarities(peptide):
    """Calculate combined similarities for CDR3a and CDR3b together."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    
    # Load full and binder data for both CDR3a and CDR3b
    cdr3a_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_full.npz"), allow_pickle=True)
    cdr3b_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_full.npz"), allow_pickle=True)
    cdr3a_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_binder.npz"), allow_pickle=True)
    cdr3b_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_binder.npz"), allow_pickle=True)
    
    results = []
    
    # Process each test sample
    for i in range(len(cdr3a_test['embeddings'])):
        test_raw_index = cdr3a_test['raw_indexes'][i]
        test_binder = cdr3a_test['binder_values'][i]
        
        best_combined_score = -2
        best_similarity_a = -1
        best_similarity_b = -1
        
        # Compare with each training sample
        for j in range(len(cdr3a_train['embeddings'])):
            train_raw_index = cdr3a_train['raw_indexes'][j]
            
            if train_raw_index == test_raw_index:
                continue
            
            # Calculate similarities for both CDR3a and CDR3b
            similarity_a = cosine_similarity(
                cdr3a_test['embeddings'][i].reshape(1, -1),
                cdr3a_train['embeddings'][j].reshape(1, -1)
            )[0][0]
            
            similarity_b = cosine_similarity(
                cdr3b_test['embeddings'][i].reshape(1, -1),
                cdr3b_train['embeddings'][j].reshape(1, -1)
            )[0][0]
            
            combined_score = similarity_a + similarity_b
            
            if combined_score > best_combined_score:
                best_combined_score = combined_score
                best_similarity_a = similarity_a
                best_similarity_b = similarity_b
        
        results.append({
            'raw_index': test_raw_index,
            'max_similarity_a': best_similarity_a,
            'max_similarity_b': best_similarity_b,
            'sum_max_similarity': best_combined_score,
            'binder': test_binder
        })
    
    # Save results
    output_file = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved combined CDR3 results to {output_file}")

# TASK 3: All CDRs Combined Analysis
#################################

def calculate_all_cdr_similarities(peptide):
    """Calculate similarities considering all CDR regions together."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    cdr_regions = ['CDR1a', 'CDR1b', 'CDR2a', 'CDR2b', 'CDR3a', 'CDR3b']
    
    # Load all test and train data
    test_data = {}
    train_data = {}
    for region in cdr_regions:
        test_data[region] = np.load(os.path.join(base_dir, f"{peptide}_{region}_full.npz"), allow_pickle=True)
        train_data[region] = np.load(os.path.join(base_dir, f"{peptide}_{region}_binder.npz"), allow_pickle=True)
    
    results = []
    
    # Process each test sample
    for i in range(len(test_data['CDR1a']['embeddings'])):
        test_raw_index = test_data['CDR1a']['raw_indexes'][i]
        test_binder = test_data['CDR1a']['binder_values'][i]
        
        best_combined_score = float('-inf')
        best_similarities = {region: -1 for region in cdr_regions}
        
        # Compare with each training sample
        for j in range(len(train_data['CDR1a']['embeddings'])):
            train_raw_index = train_data['CDR1a']['raw_indexes'][j]
            
            if train_raw_index == test_raw_index:
                continue
            
            # Calculate similarities for all regions
            current_similarities = {}
            for region in cdr_regions:
                similarity = cosine_similarity(
                    test_data[region]['embeddings'][i].reshape(1, -1),
                    train_data[region]['embeddings'][j].reshape(1, -1)
                )[0][0]
                current_similarities[region] = similarity
            
            # Calculate regular and weighted sums
            current_sum = sum(current_similarities.values())
            weighted_sum = (
                current_similarities['CDR1a'] + 
                current_similarities['CDR1b'] + 
                current_similarities['CDR2a'] + 
                current_similarities['CDR2b'] + 
                4 * current_similarities['CDR3a'] + 
                4 * current_similarities['CDR3b']
            )
            
            if weighted_sum > best_combined_score:
                best_combined_score = weighted_sum
                best_similarities = current_similarities.copy()
        
        # Store results
        result = {
            'raw_index': test_raw_index,
            'binder': test_binder,
            'weighted_sum': best_combined_score,
            'regular_sum': sum(best_similarities.values())
        }
        # Add individual similarities
        for region in cdr_regions:
            result[f'max_similarity_{region}'] = best_similarities[region]
        
        results.append(result)
    
    # Save results
    output_file = os.path.join(base_dir, f"{peptide}_all_CDR_combined_similarities.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved all CDR combined results to {output_file}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptide = sys.argv[1]
    setup_logging(peptide)
    
    try:
        # Run individual CDR3 analysis
        process_cdr3_similarities(peptide)
        logging.info("Completed individual CDR3 analysis")
        
        # Run combined CDR3 analysis
        calculate_combined_cdr3_similarities(peptide)
        logging.info("Completed combined CDR3 analysis")
        
        # Run all CDR combined analysis
        calculate_all_cdr_similarities(peptide)
        logging.info("Completed all CDR combined analysis")
        
        logging.info(f"All processing completed for peptide {peptide}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()