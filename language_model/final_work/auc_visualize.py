import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def create_combined_auc_visualization(peptides, base_path):
    """Create and save bar plot of AUC scores for multiple peptides with averages."""
    # Initialize list to store all peptides' data
    all_data = []
    
    # Mapping for the order and names we want
    name_mapping = {
        'All_CDR_Weighted': 'Weighted',
        'All_CDR_Unweighted': 'Unweighted',
        'Combined_CDR3': 'CDR3',
        'Individual_CDR3b': 'CDR3 beta',
        'Individual_CDR3a': 'CDR3 alpha'
    }
    
    # Load data for each peptide
    for peptide in peptides:
        peptide_dir = os.path.join(base_path, f"cdr_embeddings_{peptide}")
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
    
    # Reorder columns to match desired order
    desired_order = ['Weighted', 'Unweighted', 'CDR3', 'CDR3 beta', 'CDR3 alpha']
    plot_data = plot_data[desired_order]
    
    # Set up the plot
    plt.figure(figsize=(15, 8))
    
    # Create the bar plot
    x = np.arange(len(plot_data.index))
    width = 0.15
    multiplier = 0
    
    # Colors matching your example
    colors = ['#ffcdb2', '#ffb4a2', '#e5989b', '#b5838d', '#6d6875']
    
    # Plot each segment
    for attribute, color in zip(desired_order, colors):
        offset = width * multiplier
        plt.bar(x + offset, plot_data[attribute], width, label=attribute, color=color)
        multiplier += 1
    
    # Customize the plot
    plt.title('AUC Values for ESM2 embedding (sum)', fontsize=14, pad=20)
    plt.xlabel('Peptides')
    plt.ylabel('AUC Values')
    plt.ylim(0, 1.0)
    
    # Set x-axis ticks
    plt.xticks(x + width * 2, plot_data.index, rotation=0)
    
    # Add grid
    plt.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
    
    # Add legend
    plt.legend(title='Segments')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(base_path, "combined_auc_visualization.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
    base_path = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work"
    create_combined_auc_visualization(peptides, base_path)

if __name__ == "__main__":
    main()