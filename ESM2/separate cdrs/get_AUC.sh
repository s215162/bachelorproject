#!/usr/bin/env python3
import pandas as pd
from sklearn.metrics import roc_auc_score
import os

# Define directories and lists
input_directory = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim"
output_directory = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim/AUC"

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

peptides = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFKQLL",
    "YLQPRTFLL",
]

# Define analysis types and their corresponding file patterns and score columns
analysis_types = {
    "CDR3a": {
        "pattern": "similarities_{peptide}_CDR3a_all.csv",
        "score_col": "max_similarity"
    },
    "CDR3b": {
        "pattern": "similarities_{peptide}_CDR3b_all.csv",
        "score_col": "max_similarity"
    },
    "CDR3_combined": {
        "pattern": "{peptide}_CDR3_combined_all.csv",
        "score_col": "sum_max_similarity"
    },
    "Weighted": {
        "pattern": "{peptide}_all_CDR_all.csv",
        "score_col": "weighted_sum"
    },
    "Unweighted": {
        "pattern": "{peptide}_all_CDR_all.csv",
        "score_col": "unweighted_sum"
    }
}

results = []

for analysis_type, config in analysis_types.items():
    print(f"Calculating AUC for analysis type: {analysis_type}")
    
    for peptide in peptides:
        file_path = os.path.join(input_directory, config["pattern"].format(peptide=peptide))
        
        print(f"Processing file: {file_path}")
        
        try:
            # Load the data
            data = pd.read_csv(file_path)
            
            # For the all_CDR file, we can get both weighted and unweighted from same file
            if analysis_type in ["Weighted", "Unweighted"]:
                y_true = data['binder']
                y_scores = data[config["score_col"]]
                auc_score = roc_auc_score(y_true, y_scores)
                print(f"AUC score for {peptide} ({analysis_type}): {auc_score}")
                
                results.append({
                    "Peptide": peptide,
                    "Analysis_Type": analysis_type,
                    "AUC_Score": auc_score
                })
            
            # For other files, process normally
            else:
                y_true = data['binder']
                y_scores = data[config["score_col"]]
                auc_score = roc_auc_score(y_true, y_scores)
                print(f"AUC score for {peptide} ({analysis_type}): {auc_score}")
                
                results.append({
                    "Peptide": peptide,
                    "Analysis_Type": analysis_type,
                    "AUC_Score": auc_score
                })
                
        except FileNotFoundError:
            print(f"File not found for peptide {peptide}, analysis type {analysis_type}")
        except Exception as e:
            print(f"Error processing {peptide} ({analysis_type}): {e}")

# Save results
results_df = pd.DataFrame(results)

# Save detailed results
results_df.to_csv(os.path.join(output_directory, "auc_scores_detailed.csv"), index=False)

# Create pivot table for easier reading
pivot_df = results_df.pivot(index='Peptide', columns='Analysis_Type', values='AUC_Score')
pivot_df.to_csv(os.path.join(output_directory, "auc_scores_pivot.csv"))

print("AUC scores have been calculated and saved.")