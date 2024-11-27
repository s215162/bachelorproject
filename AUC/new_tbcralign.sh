#!/bin/bash

# Paths
TBCRALIGN="/home/people/morni/bin/tbcr_align"
OUTPUTDIRECTORY="/net/mimer/mnt/tank/projects2/emison/AUC/concatenated_scores"  # set output directory

# Create output directory if it doesn't exist
mkdir -p "$OUTPUTDIRECTORY"

# Define peptides
peptides=("ELAGIGILTV" "GILGFVFTL" "GLCTLVAML" "LLWNGPMAV" "RAKFKQLL" "YLQPRTFLL")

# Loop through each peptide
for peptide in "${peptides[@]}"; do
  CONCATENATED_OUTPUT="${OUTPUTDIRECTORY}/${peptide}_pred_concat_unweighted"

  # create/ clear concatenated output file
  > "$CONCATENATED_OUTPUT"

  # Loop through peptide partitions
  for i in {0..4}; do
    # Define input files
    TRAIN_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/split_data_real/train_${i}_${peptide}.csv"
    TEST_FILE="/net/mimer/mnt/tank/projects2/emison/AUC/split_data_real/test_${i}_${peptide}.csv"
    OUTPUT_FILE="${OUTPUTDIRECTORY}/test_${i}_pred_unweighted"

    # Check if the input files exist
    if [[ -f "$TRAIN_FILE" && -f "$TEST_FILE" ]]; then
      echo "Running tbcr_align for partition ${i} for peptide ${peptide}..."
      /home/people/morni/bin/tbcr_align -db $TRAIN_FILE -w "1,1,1,1,1,1" $TEST_FILE > $OUTPUT_FILE

      # Concatenate the output to the file
      cat $OUTPUT_FILE >> $CONCATENATED_OUTPUT
    else
      echo "Input files for partition $i and peptide $peptide do not exist."
    fi
  done
done
# Loop through each peptide
for peptide in "${peptides[@]}"; do

	awk '!/^#/ { print $10, $19 }' /net/mimer/mnt/tank/projects2/emison/AUC/concatenated_scores/${peptide}_pred_concat_unweighted > /net/mimer/mnt/tank/projects2/emison/AUC/concatenated_scores/${peptide}_unweighted_score_binder_only.csv

done
