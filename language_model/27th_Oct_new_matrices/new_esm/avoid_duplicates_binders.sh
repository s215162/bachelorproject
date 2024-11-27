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

# Script path for processing
script_path="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/run_esm_binders_new.py"

# Function to check if a session exists
session_exists() {
    tmux has-session -t "$1" 2>/dev/null
    return $?
}

# Loop through peptides and weights, and only create missing or duplicate-prefixed sessions
for peptide in "${peptides[@]}"; do
    for weight in "${weights[@]}"; do
        session_name="${peptide}_${weight}"
        
        # Check if the session already exists
        if session_exists "$session_name"; then
            # If session exists, create a new one prefixed with "binder_"
            new_session_name="binder_${session_name}"
            if ! session_exists "$new_session_name"; then
                # Create a new session for the duplicate case
                tmux new-session -d -s "$new_session_name" "python $script_path $peptide $weight"
                echo "Created duplicate session: $new_session_name"
            else
                echo "Session already exists: $new_session_name"
            fi
        fi
    done
done

echo "Duplicate tmux sessions created where necessary."

