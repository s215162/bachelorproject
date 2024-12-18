import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
import os
import logging
from scipy.interpolate import interp1d

def setup_logging():
    log_file = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/auc01_calculation.log"
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

def calculate_auc01(y_true, y_scores):
    """Calculate AUC0.1 score."""
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    
    interpolator = interp1d(fpr, tpr, kind='linear', bounds_error=False)
    specific_fpr = 0.1
    specific_tpr = float(interpolator(specific_fpr))
    
    mask = fpr <= specific_fpr
    auc01 = np.trapz(tpr[mask], fpr[mask]) / specific_fpr
    
    return auc01

def calculate_aucs01(peptide, base_path):
    """Calculate AUC0.1 scores for all similarity measures."""
    base_dir = os.path.join(base_path, "single_chains")
    results = []
    
    # 1. Individual CDR3 regions
    for region in ['CDR3a', 'CDR3b']:
        file_path = os.path.join(base_dir, f"{peptide}_{region}_similarities.csv")
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            try:
                auc01_score = calculate_auc01(df['binder'], df['max_similarity'])
                results.append({
                    'Peptide': peptide,
                    'Weight': f'Individual_{region}',
                    'AUC01_Score': auc01_score
                })
                logging.info(f"AUC0.1 score for {peptide} {region}: {auc01_score}")
            except Exception as e:
                logging.error(f"Error calculating AUC0.1 for {peptide} {region}: {e}")
    
    # 2. Combined CDR3 - Both sum and mean based
    # Sum-based
    cdr3_sum_path = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_sum.csv")
    if os.path.exists(cdr3_sum_path):
        df = pd.read_csv(cdr3_sum_path)
        try:
            auc01_score = calculate_auc01(df['binder'], df['sum_max_similarity'])
            results.append({
                'Peptide': peptide,
                'Weight': 'Combined_CDR3_Sum',
                'AUC01_Score': auc01_score
            })
            logging.info(f"AUC0.1 score for {peptide} Combined CDR3 (Sum): {auc01_score}")
        except Exception as e:
            logging.error(f"Error calculating AUC0.1 for {peptide} Combined CDR3 (Sum): {e}")
            
    # Mean-based
    cdr3_mean_path = os.path.join(base_dir, f"{peptide}_CDR3_combined_similarities_mean.csv")
    if os.path.exists(cdr3_mean_path):
        df = pd.read_csv(cdr3_mean_path)
        try:
            auc01_score = calculate_auc01(df['binder'], df['sum_max_similarity'])
            results.append({
                'Peptide': peptide,
                'Weight': 'Combined_CDR3_Mean',
                'AUC01_Score': auc01_score
            })
            logging.info(f"AUC0.1 score for {peptide} Combined CDR3 (Mean): {auc01_score}")
        except Exception as e:
            logging.error(f"Error calculating AUC0.1 for {peptide} Combined CDR3 (Mean): {e}")
    
    # 3. All CDRs combined
    all_cdr_path = os.path.join(base_dir, f"{peptide}_all_CDR_combined_similarities.csv")
    if os.path.exists(all_cdr_path):
        df = pd.read_csv(all_cdr_path)
        try:
            # Unweighted sum
            auc01_score_unweighted = calculate_auc01(df['binder'], df['unweighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Unweighted',
                'AUC01_Score': auc01_score_unweighted
            })
            
            # Weighted sum
            auc01_score_weighted = calculate_auc01(df['binder'], df['weighted_sum'])
            results.append({
                'Peptide': peptide,
                'Weight': 'All_CDR_Weighted',
                'AUC01_Score': auc01_score_weighted
            })
            
        except Exception as e:
            logging.error(f"Error calculating AUC0.1 for {peptide} All CDRs: {e}")
    
    # Save individual peptide results
    output_file = os.path.join(base_dir, f"{peptide}_auc01_scores.csv")
    pd.DataFrame(results).to_csv(output_file, index=False)
    logging.info(f"Saved AUC0.1 scores for {peptide} to {output_file}")
    
    return results

def create_combined_auc01_visualization(all_results, base_path):
    """Create and save bar plot of AUC0.1 scores for all peptides with averages."""
    combined_df = pd.DataFrame(all_results)
    
    name_mapping = {
        'All_CDR_Weighted': 'Weighted',
        'All_CDR_Unweighted': 'Unweighted',
        'Combined_CDR3_Sum': 'CDR3 (Sum)',
        'Combined_CDR3_Mean': 'CDR3 (Mean)',
        'Individual_CDR3b': 'CDR3 beta',
        'Individual_CDR3a': 'CDR3 alpha'
    }
    
    combined_df = combined_df[combined_df['Weight'].isin(name_mapping.keys())].copy()
    combined_df['Weight'] = combined_df['Weight'].map(name_mapping)
    
    avg_scores = combined_df.groupby('Weight')['AUC01_Score'].mean().reset_index()
    avg_scores['Peptide'] = 'AVERAGED PEPTIDE'
    
    combined_df = pd.concat([combined_df, avg_scores], ignore_index=True)
    
    plot_data = combined_df.pivot(index='Peptide', columns='Weight', values='AUC01_Score')
    
    desired_order = ['Weighted', 'Unweighted', 'CDR3 (Sum)', 'CDR3 (Mean)', 'CDR3 beta', 'CDR3 alpha']
    plot_data = plot_data[desired_order]
    
    # Set up the plot with white background
    plt.figure(figsize=(15, 8))
    plt.rcParams['figure.facecolor'] = 'white'
    
    x = np.arange(len(plot_data.index))
    width = 0.13
    multiplier = 0
    
    # Pastel colors from the image
    colors = ['#9B8BB4', '#E7B7C8', '#FFD6BA', '#E79C8E', '#A67F6C', '#B5838D']
    
    for attribute, color in zip(desired_order, colors):
        offset = width * multiplier
        plt.bar(x + offset, plot_data[attribute], width, label=attribute, color=color)
        multiplier += 1
    
    plt.title('AUC0.1 Values for ESM2 embedding', fontsize=14, pad=20)
    plt.xlabel('Peptides')
    plt.ylabel('AUC0.1 Values')
    plt.ylim(0, 1.0)
    
    plt.xticks(x + width * 2.5, plot_data.index, rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
    
    plt.legend(title='Segments', frameon=True, facecolor='white', framealpha=1)
    
    plt.tight_layout()
    
    output_path = os.path.join(base_path, "combined_auc01_visualization.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work"
    setup_logging()
    
    all_results = []
    
    try:
        for peptide in peptides:
            logging.info(f"Processing peptide: {peptide}")
            results = calculate_aucs01(peptide, base_path)
            all_results.extend(results)
        
        create_combined_auc01_visualization(all_results, base_path)
        logging.info("Processing completed for all peptides")
        
        combined_output = os.path.join(base_path, "combined_auc01_scores.csv")
        pd.DataFrame(all_results).to_csv(combined_output, index=False)
        logging.info("Saved combined results")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()