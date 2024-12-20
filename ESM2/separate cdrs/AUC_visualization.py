import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Update file paths and filenames to fit your directory

# Load the data from the file
file_path = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim/AUC/auc_scores_detailed.csv"
data = pd.read_csv(file_path)

# Define the peptides and segments
peptides = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFKQLL",
    "YLQPRTFLL"
]

segment_labels = ["Weighted", "Unweighted", "CDR3_combined", "CDR3a", "CDR3b"]

# Initialize the data matrix
auc_values = []

# Iterate through each peptide
for peptide in peptides:
    peptide_data = data[data["Peptide"] == peptide]
    auc_row = []
    for analysis_type in segment_labels:
        score = peptide_data.loc[peptide_data["Analysis_Type"] == analysis_type, "AUC_Score"]
        if not score.empty:
            auc_row.append(score.values[0])
        else:
            auc_row.append(None)
    auc_values.append(auc_row)

# Calculate averages across rows
averages_row = []
for i in range(len(segment_labels)):
    scores = [row[i] for row in auc_values if row[i] is not None]
    if scores:
        average_score = sum(scores) / len(scores)
        averages_row.append(average_score)
    else:
        averages_row.append(None)

# Append averages row
auc_values.append(averages_row)

# Add average element to list for plotting
peptide_names = peptides + ["AVERAGED PEPTIDE"]

# Create an array for the x positions of the bars
x = np.arange(len(peptide_names))
width = 0.15  # Width of the bars

# define colors for segments
colors = ["#f9cbb6", "#f3abc2", "#d7b0d2", "#c2badc", "#a3c7eb"]

# Create figure and axes
fig, ax = plt.subplots(figsize=(12, 6))

# Plot bars for each segment
for i in range(len(segment_labels)):
    ax.bar(
        x + i * width,
        [auc_values[j][i] for j in range(len(peptide_names))],
        width,
        label=segment_labels[i],
        color=colors[i]
    )

# Customize the plot
ax.set_ylim(0, 1)
ax.set_ylabel("AUC Score")
ax.set_title("AUC Values for Different Analysis Methods")
ax.set_xticks(x + width * (len(segment_labels) - 1) / 2)
ax.set_xticklabels(peptide_names, rotation=45, ha='right')
ax.legend(title="Analysis Method")

# Adjust layout and save
plt.tight_layout()
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/divided_sequences_final/cos_sim/AUC/visualization"
os.makedirs(output_dir, exist_ok=True)
output_file_path = os.path.join(output_dir, "auc_scores_comparison.png")
plt.savefig(output_file_path, format='png', bbox_inches='tight', dpi=300)
plt.close(fig)

print(f"Plot saved to {output_file_path}")

# Print the numerical values for verification
print("\nAUC Values:")
for i, peptide in enumerate(peptide_names):
    print(f"\n{peptide}:")
    for j, segment in enumerate(segment_labels):
        print(f"{segment}: {auc_values[i][j]:.3f if auc_values[i][j] is not None else 'N/A'}")