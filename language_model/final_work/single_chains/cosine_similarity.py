import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import logging
import sys
import os

def setup_logging(peptide):
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/new_similarity_processing_{peptide}.log"
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
    
    # Sum each embedding matrix along rows to get 1280-dimension vectors
    test_summed = [np.sum(emb, axis=0) for emb in test_embeddings]
    train_summed = [np.sum(emb, axis=0) for emb in train_embeddings]
    
    for i, test_vec in enumerate(test_summed):
        test_raw_index = test_raw_indexes[i]
        test_binder = test_binder_values[i]
        
        max_similarity = -1
        
        for j, train_vec in enumerate(train_summed):
            train_raw_index = train_raw_indexes[j]
            
            if train_raw_index == test_raw_index:
                continue
                
            similarity = cosine_similarity(
                test_vec.reshape(1, -1),
                train_vec.reshape(1, -1)
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
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
    
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
    """Calculate combined similarities for CDR3a and CDR3b together, using both sum and mean."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
    
    cdr3a_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_full.npz"), allow_pickle=True)
    cdr3b_test = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_full.npz"), allow_pickle=True)
    cdr3a_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3a_binder.npz"), allow_pickle=True)
    cdr3b_train = np.load(os.path.join(base_dir, f"{peptide}_CDR3b_binder.npz"), allow_pickle=True)
    
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
    
    for i in range(len(cdr3a_test_summed)):
        test_raw_index = cdr3a_test['raw_indexes'][i]
        test_binder = cdr3a_test['binder_values'][i]
        
        best_combined_score_sum = -2
        best_similarity_a_sum = -1
        best_similarity_b_sum = -1
        
        best_combined_score_mean = -2
        best_similarity_a_mean = -1
        best_similarity_b_mean = -1
        
        for j in range(len(cdr3a_train_summed)):
            train_raw_index = cdr3a_train['raw_indexes'][j]
            
            if train_raw_index == test_raw_index:
                continue
            
            # Sum-based similarities
            similarity_a_sum = cosine_similarity(
                cdr3a_test_summed[i].reshape(1, -1),
                cdr3a_train_summed[j].reshape(1, -1)
            )[0][0]
            
            similarity_b_sum = cosine_similarity(
                cdr3b_test_summed[i].reshape(1, -1),
                cdr3b_train_summed[j].reshape(1, -1)
            )[0][0]
            
            combined_score_sum = similarity_a_sum + similarity_b_sum
            
            # Mean-based similarities
            similarity_a_mean = cosine_similarity(
                cdr3a_test_mean[i].reshape(1, -1),
                cdr3a_train_mean[j].reshape(1, -1)
            )[0][0]
            
            similarity_b_mean = cosine_similarity(
                cdr3b_test_mean[i].reshape(1, -1),
                cdr3b_train_mean[j].reshape(1, -1)
            )[0][0]
            
            combined_score_mean = similarity_a_mean + similarity_b_mean
            
            if combined_score_sum > best_combined_score_sum:
                best_combined_score_sum = combined_score_sum
                best_similarity_a_sum = similarity_a_sum
                best_similarity_b_sum = similarity_b_sum
                
            if combined_score_mean > best_combined_score_mean:
                best_combined_score_mean = combined_score_mean
                best_similarity_a_mean = similarity_a_mean
                best_similarity_b_mean = similarity_b_mean
        
        results_sum.append({
            'raw_index': test_raw_index,
            'max_similarity_a': best_similarity_a_sum,
            'max_similarity_b': best_similarity_b_sum,
            'sum_max_similarity': best_combined_score_sum,
            'binder': test_binder
        })
        
        results_mean.append({
            'raw_index': test_raw_index,
            'max_similarity_a': best_similarity_a_mean,
            'max_similarity_b': best_similarity_b_mean,
            'sum_max_similarity': best_combined_score_mean,
            'binder': test_binder
        })

    # Save results   
    pd.DataFrame(results_sum).to_csv(
        os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_sum.csv"), 
        index=False
    )
    pd.DataFrame(results_mean).to_csv(
        os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_mean.csv"), 
        index=False
    )
    logging.info(f"Saved both sum and mean CDR3 results to separate files")


    logging.info(f"Saved combined CDR3 sum results to {base_dir}/{peptide}_CDR3_combined_similarities_sum.csv")
    logging.info(f"Saved combined CDR3 mean results to {base_dir}/{peptide}_CDR3_combined_similarities_mean.csv")

# TASK 3: All CDRs Combined Analysis
#################################

def calculate_all_cdr_similarities(peptide):
    """Calculate similarities considering all CDR regions together, with both weighted and unweighted sums."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
    cdr_regions = ['CDR1a', 'CDR1b', 'CDR2a', 'CDR2b', 'CDR3a', 'CDR3b']
    
    # Load all test and train data
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
    
    # Sum embeddings for all regions
    test_summed = {}
    train_summed = {}
    for region in cdr_regions:
        test_summed[region] = [np.sum(emb, axis=0) for emb in test_data[region]['embeddings']]
        train_summed[region] = [np.sum(emb, axis=0) for emb in train_data[region]['embeddings']]
    
    results = []
    
    # Process each test sample
    for i in range(len(test_summed['CDR1a'])):
        test_raw_index = test_data['CDR1a']['raw_indexes'][i]
        test_binder = test_data['CDR1a']['binder_values'][i]
        
        best_weighted_score = float('-inf')
        best_similarities = {region: -1 for region in cdr_regions}
        best_unweighted_sum = -1
        
        # Compare with each training sample
        for j in range(len(train_summed['CDR1a'])):
            train_raw_index = train_data['CDR1a']['raw_indexes'][j]
            
            if train_raw_index == test_raw_index:
                continue
            
            # Calculate similarities for all regions
            current_similarities = {}
            for region in cdr_regions:
                similarity = cosine_similarity(
                    test_summed[region][i].reshape(1, -1),
                    train_summed[region][j].reshape(1, -1)
                )[0][0]
                current_similarities[region] = similarity
            
            # Calculate unweighted sum
            unweighted_sum = (
                current_similarities['CDR1a'] + 
                current_similarities['CDR1b'] + 
                current_similarities['CDR2a'] + 
                current_similarities['CDR2b'] + 
                current_similarities['CDR3a'] + 
                current_similarities['CDR3b']
            )
            
            # Calculate weighted sum (4x weight for CDR3)
            weighted_sum = (
                current_similarities['CDR1a'] + 
                current_similarities['CDR1b'] + 
                current_similarities['CDR2a'] + 
                current_similarities['CDR2b'] + 
                4 * current_similarities['CDR3a'] + 
                4 * current_similarities['CDR3b']
            )
            
            if weighted_sum > best_weighted_score:
                best_weighted_score = weighted_sum
                best_similarities = current_similarities.copy()
                best_unweighted_sum = unweighted_sum
        
        # Store results
        result = {
            'raw_index': test_raw_index,
            'binder': test_binder,
            'weighted_sum': best_weighted_score,
            'unweighted_sum': best_unweighted_sum
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
        # Run all analyses
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