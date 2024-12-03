import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import os
import logging
import sys

def setup_logging(peptide):
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/auc_calculation_{peptide}.log"
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

def calculate_aucs(peptide):
    """Calculate AUC scores for all similarity measures."""
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    results = []
    
    # 1. Individual CDR3 regions
    for region in ['CDR3a', 'CDR3b']:
        file_path = os.path.join(base_dir, f"{peptide}_{region}_similarities.csv")
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            try:
                auc_score = roc_auc_score(df['binder'], df['max_similarity'])
                results.append({
                    'Peptide': peptide,
                    'Weight': f'Individual_{region}',
                    'AUC_Score': auc_score
                })
                logging.info(f"AUC score for {region}: {auc_score}")
            except Exception as e:
                logging.error(f"Error calculating AUC for {region}: {e}")
    
    # 2. Combined CDR3
    cdr3_path = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities.csv")
    if os.path.exists(cdr3_path):
        df = pd.read_csv(cdr3_path)
        try:
            auc_score = roc_auc_score(df['binder'], df['sum_max_similarity'])
            results.append({
                'Peptide': peptide,
                'Weight': 'Combined_CDR3',
                'AUC_Score': auc_score
            })
            logging.info(f"AUC score for Combined CDR3: {auc_score}")
        except Exception as e:
            logging.error(f"Error calculating AUC for Combined CDR3: {e}")
    
    # 3. All CDRs combined
    all_cdr_path = os.path.join(base_dir, f"{peptide}_all_CDR_combined_similarities.csv")
    if os.path.exists(all_cdr_path):
        df = pd.read_csv(all_cdr_path)
        try:
            # Unweighted sum
            auc_score_unweighted = roc_auc_score(df['binder'], df['unweighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Unweighted',
                'AUC_Score': auc_score_unweighted
            })
            logging.info(f"AUC score for All CDR Unweighted: {auc_score_unweighted}")
            
            # Weighted sum
            auc_score_weighted = roc_auc_score(df['binder'], df['weighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Weighted',
                'AUC_Score': auc_score_weighted
            })
            logging.info(f"AUC score for All CDR Weighted: {auc_score_weighted}")
            
            # Individual CDRs from the all_CDR file
            for region in ['CDR1a', 'CDR1b', 'CDR2a', 'CDR2b', 'CDR3a', 'CDR3b']:
                auc_score = roc_auc_score(df['binder'], df[f'max_similarity_{region}'])
                results.append({
                    'Peptide': peptide,
                    'Weight': f'From_All_{region}',
                    'AUC_Score': auc_score
                })
                logging.info(f"AUC score for {region} from all CDR file: {auc_score}")
                
        except Exception as e:
            logging.error(f"Error calculating AUC for All CDRs: {e}")
    
    # Save results
    output_file = os.path.join(base_dir, f"{peptide}_auc_scores.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved AUC scores to {output_file}")
    
    return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptide = sys.argv[1]
    setup_logging(peptide)
    
    try:
        results = calculate_aucs(peptide)
        logging.info("AUC calculation completed successfully")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()