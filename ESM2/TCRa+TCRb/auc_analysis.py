import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
import sys

def setup_logging(peptide):
    """Configure logging to track analysis progress for specific peptide"""
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/analysis_{peptide}.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    # Add console output for real-time monitoring
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def calculate_aucs(peptide, single_chains_dir):
    """Calculate AUC scores for different TCR similarity measures.
    
    Processes:
    1. Individual CDR3 regions (alpha and beta chains)
    2. Combined CDR3 scores (sum and mean based)
    3. All CDRs combined (weighted and unweighted)
    """
    results = []
    
    # Process individual CDR3 regions from both chains
    for region in ['CDR3a', 'CDR3b']:
        file_path = os.path.join(single_chains_dir, f"{peptide}_{region}_similarities.csv")
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
    
    # Process combined CDR3 scores using different methods
    # Sum-based combination
    cdr3_sum_path = os.path.join(single_chains_dir, f"{peptide}_CDR3_combined_similarities_sum.csv")
    if os.path.exists(cdr3_sum_path):
        df = pd.read_csv(cdr3_sum_path)
        try:
            auc_score = roc_auc_score(df['binder'], df['sum_max_similarity'])
            results.append({
                'Peptide': peptide,
                'Weight': 'Combined_CDR3_Sum',
                'AUC_Score': auc_score
            })
            logging.info(f"AUC score for Combined CDR3 (Sum): {auc_score}")
        except Exception as e:
            logging.error(f"Error calculating AUC for Combined CDR3 (Sum): {e}")
            
    # Mean-based combination
    cdr3_mean_path = os.path.join(single_chains_dir, f"{peptide}_CDR3_combined_similarities_mean.csv")
    if os.path.exists(cdr3_mean_path):
        df = pd.read_csv(cdr3_mean_path)
        try:
            auc_score = roc_auc_score(df['binder'], df['sum_max_similarity'])
            results.append({
                'Peptide': peptide,
                'Weight': 'Combined_CDR3_Mean',
                'AUC_Score': auc_score
            })
            logging.info(f"AUC score for Combined CDR3 (Mean): {auc_score}")
        except Exception as e:
            logging.error(f"Error calculating AUC for Combined CDR3 (Mean): {e}")
    
    # Process all CDRs combined (weighted and unweighted approaches)
    all_cdr_path = os.path.join(single_chains_dir, f"{peptide}_all_CDR_combined_similarities.csv")
    if os.path.exists(all_cdr_path):
        df = pd.read_csv(all_cdr_path)
        try:
            # Calculate unweighted sum AUC
            auc_score_unweighted = roc_auc_score(df['binder'], df['unweighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Unweighted',
                'AUC_Score': auc_score_unweighted
            })
            logging.info(f"AUC score for All CDR Unweighted: {auc_score_unweighted}")
            
            # Calculate weighted sum AUC
            auc_score_weighted = roc_auc_score(df['binder'], df['weighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Weighted',
                'AUC_Score': auc_score_weighted
            })
            logging.info(f"AUC score for All CDR Weighted: {auc_score_weighted}")
            
        except Exception as e:
            logging.error(f"Error calculating AUC for All CDRs: {e}")
    
    # Save compiled results
    output_file = os.path.join(single_chains_dir, f"{peptide}_auc_scores.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved AUC scores to {output_file}")
    
    return results

def create_combined_auc_visualization(peptides, single_chains_dir):
    """Create visualization comparing AUC scores across different methods and peptides.
    
    Generates a bar plot showing AUC scores for:
    - Individual CDR3 regions
    - Combined CDR3 scores
    - Full TCR analysis (weighted/unweighted)
    Includes average scores across all peptides.
    """
    logging.info("Starting visualization creation")
    
    # Verify data availability
    for peptide in peptides:
        auc_file = os.path.join(single_chains_dir, f"{peptide}_auc_scores.csv")
        if os.path.exists(auc_file):
            logging.info(f"Found: {auc_file}")
            df = pd.read_csv(auc_file)
            logging.info(f"Data for {peptide}:\n{df.to_string()}")
        else:
            logging.info(f"Missing: {auc_file}")    
    
    all_data = []
    
    # Define friendly names for plot labels
    name_mapping = {
        'All_CDR_Weighted': 'Weighted',
        'All_CDR_Unweighted': 'Unweighted',
        'Combined_CDR3_Sum': 'CDR3 (Sum)',
        'Combined_CDR3_Mean': 'CDR3 (Mean)',
        'Individual_CDR3b': 'CDR3 beta',
        'Individual_CDR3a': 'CDR3 alpha'
    }
    
    # Collect and process data for all peptides
    for peptide in peptides:
        auc_file = os.path.join(single_chains_dir, f"{peptide}_auc_scores.csv")
        
        if os.path.exists(auc_file):
            df = pd.read_csv(auc_file)
            df = df[df['Weight'].isin(name_mapping.keys())].copy()
            df['Weight'] = df['Weight'].map(name_mapping)
            df['Peptide'] = peptide
            all_data.append(df)
    
    # Combine all peptide data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Calculate and add average scores
    avg_scores = combined_df.groupby('Weight')['AUC_Score'].mean().reset_index()
    avg_scores['Peptide'] = 'AVERAGED PEPTIDE'
    combined_df = pd.concat([combined_df, avg_scores], ignore_index=True)
    
    # Prepare data for plotting
    plot_data = combined_df.pivot(index='Peptide', columns='Weight', values='AUC_Score')
    desired_order = ['Weighted', 'Unweighted', 'CDR3 (Sum)', 'CDR3 (Mean)', 'CDR3 beta', 'CDR3 alpha']
    plot_data = plot_data[desired_order]
    
    # Create and customize plot
    plt.figure(figsize=(15, 8))
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Set up bar positions and colors
    x = np.arange(len(plot_data.index))
    width = 0.13
    multiplier = 0
    colors = ["#f9cbb6", "#f3abc2", "#d7b0d2", "#c2badc", "#a3c7eb", '#77cbda']
    
    # Create grouped bars
    for attribute, color in zip(desired_order, colors):
        offset = width * multiplier
        plt.bar(x + offset, plot_data[attribute], width, label=attribute, color=color)
        multiplier += 1
    
    # Customize plot appearance
    plt.title('AUC Values for ESM2 embedding', fontsize=14, pad=20)
    plt.xlabel('Peptides')
    plt.ylabel('AUC Values')
    plt.ylim(0, 1.0)
    plt.xticks(x + width * 2.5, plot_data.index, rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
    plt.legend(title='Segments', frameon=True, facecolor='white', framealpha=1)
    plt.tight_layout()
    
    # Save final visualization
    output_path = os.path.join(single_chains_dir, "final_combined_auc_visualization.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    # Define target peptides and paths
    peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work"
    single_chains_dir = os.path.join(base_path, "single_chains")
    
    if sys.argv[1] in peptides:
        peptide = sys.argv[1]
        setup_logging(peptide)
        try:
            # Calculate AUC scores for current peptide
            results = calculate_aucs(peptide, single_chains_dir)
            logging.info("AUC calculation completed successfully")
            
            # Create visualization if processing last peptide
            if peptide == peptides[-1]:
                create_combined_auc_visualization(peptides, single_chains_dir)
                logging.info("Visualization created successfully")
                
        except Exception as e:
            logging.error(f"An error occurred: {e}")
            raise
    else:
        print(f"Peptide must be one of: {', '.join(peptides)}")
        sys.exit(1)

if __name__ == "__main__":
    main()