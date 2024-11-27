#!/bin/bash

# Define peptides and weights
peptides=(
    "ELAGIGILTV"
    "GILGFVFTL"
    "GLCTLVAML"
    "LLWNGPMAV"
    "RAKFKQLL"
    "YLQPRTFLL"
)

weights=(
    "CDR1a"
    "CDR1b"
    "CDR2a"
    "CDR2b"
    "CDR3a"
    "CDR3b"
    "TCRa"
    "TCRb"
)

# Define dataset paths
datasets=(
    "full_sequence"               # Change as needed
    "full_sequence_binders"       # Change as needed
)

# Iterate over each peptide
for peptide in "${peptides[@]}"; do
    # Iterate over each weight
    for weight in "${weights[@]}"; do
        # Iterate over each dataset
        for dataset in "${datasets[@]}"; do
            # Create a unique tmux session name
            session_name="${peptide}_${weight}_${dataset}"

            # Start a new tmux session and run the Python script
            tmux new-session -d -s "$session_name" \
                "python script.py $peptide $weight $dataset; read -p 'Press enter to exit...'"
            
            echo "Started session: $session_name"
        done
    done
done

