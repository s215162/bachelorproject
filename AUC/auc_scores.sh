#!/usr/bin/env python3

import pandas as pd
from sklearn.metrics import roc_auc_score
import os

# Define the output directory and the list of peptides
# This is a comment
output_directory = (
    "/net/mimer/mnt/tank/projects2/emison/AUC/output"  # Set your output directory here
)
peptides = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFKQLL",
    "YLQPRTFLL",
]
weights_identifier = "unweighted"  # Replace with your weights identifier
results = []

# Loop through each peptide and load the corresponding concatenated file
for peptide in peptides:
    concat_file = "/home/projects2/emison/AUC/concatenated_scores/${peptide}_score_binder_only.csv"

    # Load the concatenated file
    try:
        data = pd.read_csv(concat_file, header=None)

        # Assuming the last column is the target (binary) and the second last is the score
        # First column (target), second column (score)
        y_true = data.iloc[:, 0]  # True labels (binary)
        y_scores = data.iloc[:, 1]  # Predicted scores

        # Calculate the ROC AUC score
        auc_score = roc_auc_score(y_true, y_scores)

        print("AUC score for " + peptide + ": " + auc_score)

        # Append results to the list
        results.append(
            {"Peptide": peptide, "AUC Score": auc_score, "Weights": weights_identifier}
        )

    except FileNotFoundError:
        print(f"File for peptide {peptide} does not exist.")
    except Exception as e:
        print(f"An error occurred while processing {peptide}: {e}")

# Save results to a text file
results_df = pd.DataFrame(results)
results_df.to_csv(
    os.path.join(output_directory, "auc_scores.txt"), index=False, sep="\t"
)

print("AUC scores have been calculated and saved.")
