import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import logging
import sys
import os
# update all paths to work for your directory

def setup_logging(peptide):
    """Configure logging for similarity analysis of TCR sequences"""
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim"
    log_file = f"{base_path}/logs/similarity_processing_{peptide}.log"
    
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

def load_embeddings(file_path):
    """Load ESM embeddings and their corresponding indices from CSV"""
    if not os.path.exists(file_path):
        logging.error(f"File not found: {file_path}")
        raise FileNotFoundError(f"Could not find file: {file_path}")
    
    df = pd.read_csv(file_path)
    raw_indexes = df['raw_index'].values
    
    embedding_columns = [col for col in df.columns if col != 'raw_index']
    embeddings = df[embedding_columns].values
    
    return embeddings, raw_indexes

def get_binder_info(raw_indices, current_peptide):
    """Map raw indices to their corresponding binder values"""
    full_data = pd.read_csv("/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv")
    full_data = full_data[full_data['peptide_x'] == current_peptide]
    binder_map = dict(zip(full_data['raw_index'], full_data['binder']))
    return [binder_map.get(idx) for idx in raw_indices]

def calculate_similarity_scores(test_emb, train_emb):
    """Calculate cosine similarity between two ESM embedding vectors"""
    return cosine_similarity(test_emb.reshape(1, -1), train_emb.reshape(1, -1))[0][0]

def process_single_region(peptide, region):
    """Process similarities for one CDR region"""
    logging.info(f"Processing {region} for peptide {peptide}")
    
    base_path = f"/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/{region.lower()}"
    output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim"
    os.makedirs(output_dir, exist_ok=True)
    
    # load test and trainibg data
    test_file = f"{base_path}/sequence_summed_vectors_swaps_{peptide}_{region}_np.csv"
    train_file = f"{base_path}/sequence_summed_vectors_binders_{peptide}_{region}_np.csv"
    
    test_emb, test_idx = load_embeddings(test_file)
    train_emb, train_idx = load_embeddings(train_file)
    
    # Get binding information
    test_binder = get_binder_info(test_idx, peptide)
    train_binder = get_binder_info(train_idx, peptide)
    
    # Get similarities
    results = []
    for i, (test_embedding, test_raw_idx, test_is_binder) in enumerate(zip(test_emb, test_idx, test_binder)):
        max_similarity = -1
        for j, train_embedding in enumerate(train_emb):
            similarity = calculate_similarity_scores(test_embedding, train_embedding)
            max_similarity = max(max_similarity, similarity)
            
        results.append({
            'raw_index': test_raw_idx,
            'max_similarity': max_similarity,
            'binder': test_is_binder
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    binders_df = results_df[results_df['binder'] == 1]
    
    binders_df.to_csv(f"{output_dir}/similarities_{peptide}_{region}_binders.csv", index=False)
    results_df.to_csv(f"{output_dir}/similarities_{peptide}_{region}_all.csv", index=False)
    
    logging.info(f"Saved results to {output_dir}")

def process_combined_cdr3(peptide):
    """Process combined CDR3a and CDR3b similarities"""
    logging.info(f"Processing combined CDR3 for peptide {peptide}")
    
    regions = ['CDR3a', 'CDR3b']
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11"
    output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load all data
    data = {}
    for region in regions:
        test_file = f"{base_path}/{region.lower()}/sequence_summed_vectors_swaps_{peptide}_{region}_np.csv"
        train_file = f"{base_path}/{region.lower()}/sequence_summed_vectors_binders_{peptide}_{region}_np.csv"
        
        data[region] = {
            'test': load_embeddings(test_file),
            'train': load_embeddings(train_file)
        }
    
    # Process similarities
    results = []
    test_embeddings_a, test_idx = data['CDR3a']['test']
    test_binder = get_binder_info(test_idx, peptide)
    
    for i, (test_idx_i, binder) in enumerate(zip(test_idx, test_binder)):
        best_combined_score = -2
        best_scores = {'CDR3a': -1, 'CDR3b': -1}
        
        for j in range(len(data['CDR3a']['train'][0])):
            current_scores = {}
            for region in regions:
                test_emb = data[region]['test'][0][i]
                train_emb = data[region]['train'][0][j]
                current_scores[region] = calculate_similarity_scores(test_emb, train_emb)
            
            combined_score = sum(current_scores.values())
            if combined_score > best_combined_score:
                best_combined_score = combined_score
                best_scores = current_scores.copy()
        
        results.append({
            'raw_index': test_idx_i,
            'max_similarity_CDR3a': best_scores['CDR3a'],
            'max_similarity_CDR3b': best_scores['CDR3b'],
            'sum_max_similarity': best_combined_score,
            'binder': binder
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    binders_df = results_df[results_df['binder'] == 1]
    
    binders_df.to_csv(f"{output_dir}/{peptide}_CDR3_combined_binders.csv", index=False)
    results_df.to_csv(f"{output_dir}/{peptide}_CDR3_combined_all.csv", index=False)
    
    logging.info(f"Saved combined CDR3 results to {output_dir}")

def process_all_regions(peptide):
    """Process all CDR regions with weighted combinations"""
    logging.info(f"Processing all CDR regions for peptide {peptide}")
    
    regions = ['CDR1a', 'CDR1b', 'CDR2a', 'CDR2b', 'CDR3a', 'CDR3b']
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11"
    output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data for all regions
    data = {}
    for region in regions:
        test_file = f"{base_path}/{region.lower()}/sequence_summed_vectors_swaps_{peptide}_{region}_np.csv"
        train_file = f"{base_path}/{region.lower()}/sequence_summed_vectors_binders_{peptide}_{region}_np.csv"
        
        data[region] = {
            'test': load_embeddings(test_file),
            'train': load_embeddings(train_file)
        }
    
    # Define weights for different regions
    weights = {region: 4.0 if region.startswith('CDR3') else 1.0 for region in regions}
    
    # Process similarities
    results = []
    test_embeddings, test_idx = data['CDR1a']['test']
    test_binder = get_binder_info(test_idx, peptide)
    
    for i, (test_idx_i, binder) in enumerate(zip(test_idx, test_binder)):
        best_weighted_score = float('-inf')
        best_scores = {region: -1 for region in regions}
        best_unweighted_sum = -1
        
        for j in range(len(data['CDR1a']['train'][0])):
            current_scores = {}
            for region in regions:
                test_emb = data[region]['test'][0][i]
                train_emb = data[region]['train'][0][j]
                current_scores[region] = calculate_similarity_scores(test_emb, train_emb)
            
            weighted_sum = sum(score * weights[region] for region, score in current_scores.items())
            unweighted_sum = sum(current_scores.values())
            
            if weighted_sum > best_weighted_score:
                best_weighted_score = weighted_sum
                best_scores = current_scores.copy()
                best_unweighted_sum = unweighted_sum
        
        result = {
            'raw_index': test_idx_i,
            'weighted_sum': best_weighted_score,
            'unweighted_sum': best_unweighted_sum,
            'binder': binder
        }
        result.update({f'max_similarity_{region}': best_scores[region] for region in regions})
        results.append(result)
    
    # Save results
    results_df = pd.DataFrame(results)
    binders_df = results_df[results_df['binder'] == 1]
    
    binders_df.to_csv(f"{output_dir}/{peptide}_all_CDR_binders.csv", index=False)
    results_df.to_csv(f"{output_dir}/{peptide}_all_CDR_all.csv", index=False)
    
    logging.info(f"Saved all CDR results to {output_dir}")

def main(peptide):
    """Main function to process all similarity calculations"""
    setup_logging(peptide)
    
    try:
        # Process each analysis type
        for region in ['CDR3a', 'CDR3b']:
            process_single_region(peptide, region)
        
        process_combined_cdr3(peptide)
        process_all_regions(peptide)
        
        logging.info(f"All processing completed for peptide {peptide}")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    if len(sys.argv) != 2:
        logging.error("Incorrect number of arguments")
        sys.exit(1)
    
    peptide = sys.argv[1]
    main(peptide)