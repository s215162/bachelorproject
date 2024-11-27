#!/bin/bash

# Paths
TBCRALIGN="/home/people/morni/bin/tbcr_align"
OUTPUTDIRECTORY="./output_testing"  # Set your output directory here

# Create output directory if it doesn't exist
mkdir -p "$OUTPUTDIRECTORY"

TRAIN_FILE="/net/mimer/mnt/tank/projects2/emison/dist_matrix1/train.csv"
TEST_FILE="/net/mimer/mnt/tank/projects2/emison/dist_matrix1/test.csv"
OUTPUT_FILE="${OUTPUTDIRECTORY}/test_pred"

$TBCRALIGN -db "$TRAIN_FILE" "$TEST_FILE" > "$OUTPUT_FILE"




