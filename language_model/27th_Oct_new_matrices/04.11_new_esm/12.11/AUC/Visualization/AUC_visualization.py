import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Load the data from the file
file_path = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/redo/auc_scores.csv"
data = pd.read_csv(file_path, sep=",")  

# Strip spaces from column names
data.columns = data.columns.str.strip()

# Check if 'Peptide' column exists
try:
    unique_peptides = data["Peptide"].unique()
    print("Unique peptides:", unique_peptides)
except KeyError as e:
    print(f"KeyError: {e}")
    print("Available columns:", data.columns)

# Initialize the data matrix
auc_values = []

# Iterate through each unique peptide
for peptide in unique_peptides:
    # Filter the data for the current peptide
    peptide_data = data[data["Peptide"] == peptide]

    # Extract the AUC scores for each weight type
    scores = peptide_data["AUC Score"].tolist()

    # Ensure the order is correct (weighted, unweighted, CDR3, CDR3_A, CDR3_B)
    weights = ["Weighed", "Unweighed", "CDR3", "CDR2", "CDR1", "full"]  # Ensure weights match case
    auc_row = []

    for weight in weights:
        score = peptide_data.loc[peptide_data["Weights"] == weight, "AUC Score"]
        if not score.empty:
            auc_row.append(score.values[0])  # Get the first score if it exists
        else:
            auc_row.append(None)  # If no score exists for this weight

    auc_values.append(auc_row)  # Add the row to the data matrix

# Calculate the averages across the rows
averages_row = []
for i in range(len(weights)):
    # Extract the scores for this weight across all peptides
    scores = [row[i] for row in auc_values if row[i] is not None]

    if scores:  # Check if there are any scores to average
        average_score = sum(scores) / len(scores)
        averages_row.append(average_score)
    else:
        averages_row.append(None)  # If no scores are available, append None

# Append the averages row to auc_values
auc_values.append(averages_row)

# Display the resulting data matrix
print("AUC values matrix:")
print(auc_values)

# Labels for the segments
segment_labels = ["Weighted", "Unweighted", "CDR3", "CDR2", "CDR1", "TCRa + TCRb"]

# Names for each peptide
peptide_names = [
    "ELAGIGILTV",
    "GILGFVFTL",
    "GLCTLVAML",
    "LLWNGPMAV",
    "RAKFQLL",
    "YLQPRTFLL",
    "AVERAGED PEPTIDE",
]

# Create an array for the x positions of the bars
num_segments = len(segment_labels)
x = np.arange(len(peptide_names))  # The peptide label locations
width = 0.15  # Width of the bars

# Custom colors for each segment
ly_colors = ["#f9cbb6", "#f3abc2", "#d7b0d2", "#c2badc", "#a3c7eb", "#ddd5f3"]
mots_colors = ["#f2b22e", "#307bbc", "#197275", "#2e9adc", "#f29930"]
yofo_colors = ["#0d2680", "#f3dc1e", "#ea892d", "#aacf34", "#11c2ee"]

# Create a figure and axes
fig, ax = plt.subplots(figsize=(12, 6))

for i in range(num_segments):
    # Shift the bars to the right for each segment
    ax.bar(
        x + i * width,
        [auc_values[j][i] for j in range(len(peptide_names))],
        width,
        label=segment_labels[i],
        color=ly_colors[i],
    )

# Set y-axis limit
ax.set_ylim(0, 1)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel("AUC 0.1 for ESM2 mean ")
ax.set_title("AUC 0.1 Values for Each Peptide Segment")
ax.set_xticks(x + width * (num_segments - 1) / 2)
ax.set_xticklabels(peptide_names)
ax.legend(title="Segments")

# Save the plot as a PNG file
output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/Visualization/"  # Change to your desired output directory
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist
output_file_path = os.path.join(output_dir, "auc_scores_ESM_sum.png")
plt.tight_layout()
plt.savefig(output_file_path, format='png')
plt.close(fig)  # Close the figure to free up memory
print(f"Plot saved to {output_file_path}")

