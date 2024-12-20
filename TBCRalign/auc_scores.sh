import pandas as pd
from sklearn.metrics import roc_auc_score
import os
import sys

def calculate_auc_scores(peptides, weights_identifier, paths):
   """Calculate AUC scores for different peptides"""
   results = []
   
   for peptide in peptides:
       input_file = os.path.join(paths['input_dir'], f"{peptide}_score_binder_only.csv")
       
       try:
           data = pd.read_csv(input_file, header=None)
           y_true = data.iloc[:, 0]
           y_scores = data.iloc[:, 1]
           
           auc_score = roc_auc_score(y_true, y_scores)
           print(f"AUC score for {peptide}: {auc_score}")
           
           results.append({
               "Peptide": peptide,
               "AUC Score": auc_score,
               "Weights": weights_identifier
           })
           
       except FileNotFoundError:
           print(f"File for peptide {peptide} does not exist.")
       except Exception as e:
           print(f"An error occurred while processing {peptide}: {e}")
   
   return results

def main():
   # Add your paths here
   paths = {
       'input_dir': '',  # Directory with score files
       'output_dir': ''  # Where to save results
   }
   
   peptides = [
       "ELAGIGILTV", "GILGFVFTL", "GLCTLVAML",
       "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"
   ]
   
   weights_identifier = "unweighted"
   
   try:
       os.makedirs(paths['output_dir'], exist_ok=True)
       results = calculate_auc_scores(peptides, weights_identifier, paths)
       
       # Save results
       results_df = pd.DataFrame(results)
       output_file = os.path.join(paths['output_dir'], "auc_scores.txt")
       results_df.to_csv(output_file, index=False, sep="\t")
       print("AUC scores calculated and saved.")
       
   except Exception as e:
       print(f"An error occurred: {e}")
       sys.exit(1)

if __name__ == "__main__":
   main()