#!/bin/bash

# Path to input data file
INPUT_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/data/240905_nettcr_exp_paired_easy7peps_withswaps.txt"

# Directory to store the output files
OUTPUT_DIR="/net/mimer/mnt/tank/projects2/emison/AUC/split_data_real"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Define the TCR chains
string_chains="['A1', 'A2', 'A3', 'B1', 'B2', 'B3']"

# Loop through each unique peptide in the input data
for peptide in $(cut -d',' -f1 "$INPUT_FILE" | tail -n +2 | sort | uniq); do
    echo "Processing peptide: $peptide"

    # Loop over partitions (0 to 4)
    for partition in {0..4}; do
        echo "Processing partition: $partition for peptide: $peptide"

        # Python code to save files in the correct format
        /usr/bin/python3 -W ignore <<EOF
import os
import pandas as pd

# Define columns for tbcr_align format
chains = ${string_chains}  # Directly use the Python list format

# Load the input data
df = pd.read_csv("$INPUT_FILE")

# Filter data for the current peptide
df_peptide = df[df['peptide'] == '${peptide}']  # Pass 'peptide' from Bash

# Split the data into test and train sets based on partition
partition = ${partition}  # Use 'partition' from the Bash loop
test_df = df_peptide[df_peptide['partition'] == partition]
train_df = df_peptide[df_peptide['partition'] != partition]

# Only keep binders for training data
train_df = train_df[train_df['binder'] == 1]

# Check if the 'binder' or 'target' column exists; if neither exists, create a fake 'binder' column
binder_column = 'binder' if 'binder' in df.columns else 'target' if 'target' in df.columns else None
if binder_column is None:
    binder_column = 'binder'
    df[binder_column] = 0.5  # Create a fake column with a default value

# Select only the relevant chain columns and the binder/target column
test_df = test_df[chains + [binder_column]]
train_df = train_df[chains + [binder_column]]

# Define the output directory and file names for the test and train sets
output_dir = "$OUTPUT_DIR"
test_file = f"{output_dir}/test_{partition}_${peptide}.csv"  # Format the file name with 'peptide' and 'partition'
train_file = f"{output_dir}/train_{partition}_${peptide}.csv"

# Save the test and train sets as CSVs with tab-separated values, no headers, and no indices
test_df.to_csv(test_file, sep='\t', index=False, header=False)
train_df.to_csv(train_file, sep='\t', index=False, header=False)

print(f"Files saved for partition {partition} and peptide ${peptide}")
EOF

        # Add a new column at the beginning of each file with 'XX'
        gawk '{print "XX",$0}' "$OUTPUT_DIR/test_${partition}_${peptide}.csv" > "$OUTPUT_DIR/test_${partition}_${peptide}_updated.csv"
        gawk '{print "XX",$0}' "$OUTPUT_DIR/train_${partition}_${peptide}.csv" > "$OUTPUT_DIR/train_${partition}_${peptide}_updated.csv"

        # Remove the old files and rename the updated ones
        mv "$OUTPUT_DIR/test_${partition}_${peptide}_updated.csv" "$OUTPUT_DIR/test_${partition}_${peptide}.csv"
        mv "$OUTPUT_DIR/train_${partition}_${peptide}_updated.csv" "$OUTPUT_DIR/train_${partition}_${peptide}.csv"

        echo "Finished processing partition: $partition for peptide: $peptide"
        
    done

    echo "All partitions processed for peptide: $peptide"
done

