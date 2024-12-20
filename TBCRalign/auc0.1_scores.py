import pandas as pd
from sklearn.metrics import roc_auc_score
import os
import sys

def calculate_auc_scores(peptides, weights, paths):
    """Calcualte AUC scores for diffrent weight configurations"""
    results = []
    
    for weight in weights:
        print(f"Calculating AUC0.1 for peptides with weights: {weight}")
        
        for peptide in peptides:
            input_file = os.path.join(paths['input_dir'], f"{peptide}_{weight}_score_binder_only.csv")
            
            try:
                data = pd.read_csv(input_file, header=None, sep="\s+")
                y_true = data.iloc[:, 0]
                y_scores = data.iloc[:, 1]
                
                auc_score = roc_auc_score(y_true, y_scores, max_fpr=0.1)
                print(f"AUC0.1 score for {peptide} for {weight}: {auc_score}")
                
                results.append({
                    "Peptide": peptide,
                    "AUC0.1 Score": auc_score,
                    "Weights": weight,
                })
                
            except FileNotFoundError:
                print(f"File for peptide {peptide} does not exist.")
            except Exception as e:
                print(f"An error occurred while processing {peptide}: {e}")
    
    return results

def main():
    # Add your input/output directories here
    paths = {
        'input_dir': '',  # Directory with score files
        'output_dir': ''  # Where to save results
    }
    
    peptides = [
        "ELAGIGILTV", "GILGFVFTL", "GLCTLVAML",
        "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"
    ]
    
    weights = [
        "weighted", "unweighted", "CDR3",
        "CDR3_A", "CDR3_B"
    ]
    
    try:
        os.makedirs(paths['output_dir'], exist_ok=True)
        results = calculate_auc_scores(peptides, weights, paths)
        
        # Save results
        results_df = pd.DataFrame(results)
        output_file = os.path.join(paths['output_dir'], "auc0.1_scores.txt")
        results_df.to_csv(output_file, index=False, sep="\t")
        print("AUC scores have been calculated and saved.")
        
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()