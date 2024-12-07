import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def create_model_comparison_plot(base_path):
    """Create bar plot comparing weighted AUC scores between ESM2 and tbcralign"""
    # Load tBCR data
    tbcr_path = os.path.join(base_path, "tbcr_auc.txt")
    tbcr_df = pd.read_csv(tbcr_path)
    
    # Load ESM2 data
    esm2_results = []
    for peptide in ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]:
        auc_file = os.path.join(base_path, f"{peptide}_auc_scores.csv")
        if os.path.exists(auc_file):
            df = pd.read_csv(auc_file)
            weighted_row = df[df['Weight'] == 'All_CDR_Weighted'].iloc[0]
            esm2_results.append({
                'Peptide': peptide,
                'AUC_Score': weighted_row['AUC_Score'],
                'Model': 'ESM2'
            })
    
    esm2_df = pd.DataFrame(esm2_results)
    # Add column to tbcr data for easier distinction between the two models
    tbcr_df['Model'] = 'tBCR'
    
    # Calculate averages
    esm2_avg = pd.DataFrame([{
        'Peptide': 'AVERAGED PEPTIDE',
        'AUC_Score': esm2_df['AUC_Score'].mean(),
        'Model': 'ESM2'
    }])
    
    tbcr_avg = pd.DataFrame([{
        'Peptide': 'AVERAGED PEPTIDE',
        'AUC_Score': tbcr_df['AUC_Score'].mean(),
        'Model': 'tBCR'
    }])
    
    # Combine all data
    plot_df = pd.concat([esm2_df, tbcr_df, esm2_avg, tbcr_avg], ignore_index=True)
    
    # Create the plot
    plt.figure(figsize=(12, 6))
    plt.rcParams['figure.facecolor'] = 'white'
    
    # Set up positions for bars
    peptides = plot_df['Peptide'].unique()
    x = np.arange(len(peptides))
    width = 0.35
    
    # Plot bars for each model
    plt.bar(x - width/2, 
            plot_df[plot_df['Model'] == 'ESM2']['AUC_Score'],
            width, 
            label='ESM2',
            color='#9B8BB4')
    
    plt.bar(x + width/2, 
            plot_df[plot_df['Model'] == 'tBCR']['AUC_Score'],
            width, 
            label='tBCR',
            color='#E7B7C8')
    
    plt.title('Model Comparison: Weighted AUC Scores', fontsize=14, pad=20)
    plt.xlabel('Peptides')
    plt.ylabel('AUC Score')
    plt.ylim(0, 1.0)
    
    plt.xticks(x, peptides, rotation=0)
    plt.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
    
    plt.legend(frameon=True, facecolor='white', framealpha=1)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(base_path, "model_comparison.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
    create_model_comparison_plot(base_path)

if __name__ == "__main__":
    main()