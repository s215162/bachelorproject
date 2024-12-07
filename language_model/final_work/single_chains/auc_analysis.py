import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
import sys

def setup_logging(peptide):
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/analysis_{peptide}.log"
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
    base_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
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
    
    # 2. Combined CDR3 - Both sum and mean based
    # Sum-based
    cdr3_sum_path = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_sum.csv")
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
            
    # Mean-based
    cdr3_mean_path = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_mean.csv")
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
            
        except Exception as e:
            logging.error(f"Error calculating AUC for All CDRs: {e}")
    
    # Save results
    output_file = os.path.join(base_dir, f"{peptide}_auc_scores.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved AUC scores to {output_file}")
    
    return results

def create_combined_auc_visualization(peptides, base_path):
    """Create and save bar plot of AUC scores for multiple peptides with averages."""
    all_data = []
    
    name_mapping = {
        'All_CDR_Weighted': 'Weighted',
        'All_CDR_Unweighted': 'Unweighted',
        'Combined_CDR3_Sum': 'CDR3 (Sum)',
        'Combined_CDR3_Mean': 'CDR3 (Mean)',
        'Individual_CDR3b': 'CDR3 beta',
        'Individual_CDR3a': 'CDR3 alpha'
    }
    
    # Load data for each peptide
    for peptide in peptides:
        peptide_dir = os.path.join(base_path, "single_chains")
        auc_file = os.path.join(peptide_dir, f"{peptide}_auc_scores.csv")
        
        if os.path.exists(auc_file):
            df = pd.read_csv(auc_file)
            df = df[df['Weight'].isin(name_mapping.keys())].copy()
            df['Weight'] = df['Weight'].map(name_mapping)
            df['Peptide'] = peptide
            all_data.append(df)
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Calculate average scores across peptides
    avg_scores = combined_df.groupby('Weight')['AUC_Score'].mean().reset_index()
    avg_scores['Peptide'] = 'AVERAGED PEPTIDE'
    
    # Add averages to the combined dataframe
    combined_df = pd.concat([combined_df, avg_scores], ignore_index=True)
    
    # Create pivot table for plotting
    plot_data = combined_df.pivot(index='Peptide', columns='Weight', values='AUC_Score')
    
    # Define column order
    desired_order = ['Weighted', 'Unweighted', 'CDR3 (Sum)', 'CDR3 (Mean)', 'CDR3 beta', 'CDR3 alpha']
    plot_data = plot_data[desired_order]
    
    # Set up the plot with white background
    plt.figure(figsize=(15, 8))
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Create the bar plot
    x = np.arange(len(plot_data.index))
    width = 0.13
    multiplier = 0
    
    # Pastel colors from the image
    colors = ['#9B8BB4', '#E7B7C8', '#FFD6BA', '#E79C8E', '#A67F6C', '#B5838D']
    
    # Plot each segment
    for attribute, color in zip(desired_order, colors):
        offset = width * multiplier
        plt.bar(x + offset, plot_data[attribute], width, label=attribute, color=color)
        multiplier += 1
    
    plt.title('AUC Values for ESM2 embedding', fontsize=14, pad=20)
    plt.xlabel('Peptides')
    plt.ylabel('AUC Values')
    plt.ylim(0, 1.0)
    
    plt.xticks(x + width * 2.5, plot_data.index, rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
    
    # Customize legend
    plt.legend(title='Segments', frameon=True, facecolor='white', framealpha=1)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(base_path, "combined_auc_visualization.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work"
    
    if sys.argv[1] in peptides:
        peptide = sys.argv[1]
        setup_logging(peptide)
        try:
            results = calculate_aucs(peptide)
            logging.info("AUC calculation completed successfully")
            
            # If this is the last peptide, create visualization
            if peptide == peptides[-1]:
                create_combined_auc_visualization(peptides, base_path)
                logging.info("Visualization created successfully")
                
        except Exception as e:
            logging.error(f"An error occurred: {e}")
            raise
    else:
        print(f"Peptide must be one of: {', '.join(peptides)}")
        sys.exit(1)

if __name__ == "__main__":
    main()