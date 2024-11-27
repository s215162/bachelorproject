import pandas as pd
from sklearn.metrics import roc_auc_score
import os

# Define the output directory and the list of peptides
output_directory = (
    "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/Calculated_AUC"
)
peptides = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFKQLL",
    "YLQPRTFLL",
]
weights_identifier = [
    "Weighed",
    "Unweighed",
    "CDR3",
    "CDR2",
    "CDR1",
    "full"
]

results = []

for weight in weights_identifier:
    print(f"Calculating AUC for peptides with weights: {weight}")
    for peptide in peptides:
        concat_file = os.path.join(
            "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/For_AUC/databehandling/",
            f"{peptide}_{weight}.csv"
        )
        print(f"Looking for file: {concat_file}")
        try:
            # Load the concatenated file 
            data = pd.read_csv(concat_file)

            # Check columns to ensure correct indexing
            print(f"Columns in data: {data.columns.tolist()}")  # Debugging output

            # Ensure you adjust these indices based on your actual file structure
            y_true = data['binder']  # Binary true/false
            y_scores = data['sum_max_similarity']  # Similarity data

            # Calculate the ROC AUC score
            auc_score = roc_auc_score(y_true, y_scores, max_fpr=0.1)

            print(f"AUC score for {peptide} for {weight}: {auc_score}")

            results.append({
                "Peptide": peptide,
                "AUC Score": auc_score,
                "Weights": weight,
            })

        except FileNotFoundError:
            print(f"File for peptide {peptide} does not exist.")
        except Exception as e:
            print(f"An error occurred while processing {peptide}: {e}")

# Save results to a text file
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_directory, "auc0.1_scores.csv"), index=False)

print("AUC0.1 scores have been calculated and saved.")


