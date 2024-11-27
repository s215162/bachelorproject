#!/usr/bin/env python3

import pandas as pd
import sys

# Check if the correct number of arguments is provided
if len(sys.argv) != 4:
    print(
        "Usage: ./merge_sequences.py <original_data.csv> <full_sequence_data.txt> <output_file.csv>"
    )
    sys.exit(1)

# File paths from command-line arguments
original_data_file = sys.argv[1]
full_sequence_data_file = sys.argv[2]
output_file = sys.argv[3]

# Load old data
old_data = pd.read_csv(original_data_file)

# Load new data with full sequences
new_data = pd.read_csv(
    full_sequence_data_file,
    sep=" ",
    header=None,
    names=[
        "raw_index",
        "TCRa",
        "TCRb",
        "peptide",
        "CDR1a",
        "CDR2a",
        "CDR3a",
        "CDR1b",
        "CDR2b",
        "CDR3b",
    ],
)

# Merge the datasets on the 'raw_index' column
merged_data = pd.merge(old_data, new_data, on="raw_index", how="left")

# Add a new column 'tcr_full' that concatenates TCRa and TCRb
merged_data["tcr_full"] = merged_data["TCRa"] + merged_data["TCRb"]

# Save the merged data to a new file
merged_data.to_csv(output_file, index=False)

print(f"Merged data saved to {output_file}")

