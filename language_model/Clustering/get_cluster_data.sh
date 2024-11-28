#!/bin/bash

# Directory containing the input files
INPUT_DIR="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/AUC/redo"  # Change this to your files' directory
OUTPUT_DIR="/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity"
OUTPUT_FILE="$OUTPUT_DIR/tcra_tcrb_similarity.csv"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize the output file and add the header from the first file
first_file=$(ls $INPUT_DIR/*_full.csv | head -n 1)
head -n 1 "$first_file" > "$OUTPUT_FILE"

# Append all files, skipping their headers
for file in $INPUT_DIR/*_full.csv; do
    tail -n +2 "$file" >> "$OUTPUT_FILE"
done

echo "Combined file saved to: $OUTPUT_FILE"
