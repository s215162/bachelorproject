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

# Script path
script_path="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/run_esm_new.py"

# Iterate over each peptide and weight to create a tmux session
for peptide in "${peptides[@]}"; do
    for weight in "${weights[@]}"; do
        # Create a unique session name for each peptide-weight combination
        session_name="${peptide}_${weight}"
        
        # Create a new tmux session, run the Python script, and detach
        tmux new-session -d -s "$session_name" "python $script_path $peptide $weight"
        
        echo "Created tmux session: $session_name"
    done
done

echo "All tmux sessions created."

