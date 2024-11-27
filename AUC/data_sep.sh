#!/bin/bash

# Path to input data file
INPUT_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/data/240905_nettcr_exp_paired_easy7peps_withswaps.txt"

# Directory to store the output files
OUTPUT_DIR="/net/mimer/mnt/tank/projects2/emison/AUC/split_data"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each unique peptide in the input data
for peptide in $(cut -d',' -f1 "$INPUT_FILE" | tail -n +2 | sort | uniq); do
    echo "Processing peptide: $peptide"

    # Loop over the 5 partitions (0 to 4)
    for i in {0..4}; do
        echo "Processing partition: $i for peptide: $peptide"

        # Use Python to filter the data and create test/train files
        python3 <<EOF
import os
import pandas as pd

# Read input data
df = pd.read_csv("$INPUT_FILE")

# Filter data for the current peptide
df_peptide = df[df['peptide'] == '$peptide']

# Split the data into test and train sets
partition = $i  # Pass 'i' from the Bash loop
test_df = df_peptide[df_peptide['partition'] == partition]
train_df = df_peptide[df_peptide['partition'] != partition]

# Only include binders for training
train_df = train_df[train_df['binder'] == 1]

# Define the columns required for tbcr_align
columns = ['raw_index', 'A1', 'A2', 'A3', 'B1', 'B2', 'B3', 'binder']

# Select and reorder columns for tbcr_align format
test_df = test_df[columns]
train_df = train_df[columns]

# Get the output directory from the environment
output_dir = "$OUTPUT_DIR"  # Pass 'OUTPUT_DIR' from Bash

# Define file paths for test and train sets
test_file = f"{output_dir}/test_{partition}_$peptide.csv"
train_file = f"{output_dir}/train_{partition}_$peptide.csv"

# Save the test and train sets as CSVs without headers or indices
test_df.to_csv(test_file, index=False, header=False)
train_df.to_csv(train_file, index=False, header=False)
EOF

        echo "Finished processing partition: $i for peptide: $peptide"
    done

    echo "All partitions processed for peptide: $peptide"
done

