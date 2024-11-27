#!/bin/bash

# Paths
TBCRALIGN="/home/people/morni/bin/tbcr_align"
OUTPUTDIRECTORY="./output_final"  # Set your output directory here

# Create output directory if it doesn't exist
mkdir -p "$OUTPUTDIRECTORY"

# Loop through your partitions
for i in {0..4}; do
  # Define your input files
  TRAIN_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/split_data_real/train_${i}_ELAGIGILTV.csv"
  TEST_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/split_data_real/test_${i}_ELAGIGILTV.csv"
  OUTPUT_FILE="${OUTPUTDIRECTORY}/test_${i}_pred"

  # Check if the input files exist
  if [[ -f "$TRAIN_FILE" && -f "$TEST_FILE" ]]; then
    # Run tbcralign without specifying the BLOSUM file
    echo "Running tbcr_align for partition $i..."
    $TBCRALIGN -db $TRAIN_FILE $TEST_FILE > $OUTPUT_FILE
  else
    echo "Input files for partition $i do not exist."
  fi
done
