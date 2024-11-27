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

# Set the directory where the batch scripts are located
batch_script_dir="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/batch_scripts"

# Script path
script_path="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/run_esm_new.py"

# Ensure the batch_scripts directory exists
mkdir -p "$batch_script_dir"

# Iterate over each peptide and weight to create a tmux session script
for peptide in "${peptides[@]}"; do
    for weight in "${weights[@]}"; do
        # Create a unique session name for each peptide-weight combination
        session_name="${peptide}_${weight}"
        
        # Create the batch script file for the session
        batch_script="${batch_script_dir}/${session_name}_tmux.sh"
        
        # Write the tmux commands to the batch script file
        cat <<EOF > "$batch_script"
#!/bin/bash
# This script runs the Python script in a tmux session

# Start a new tmux session and run the Python script
tmux new-session -d -s "$session_name" "python $script_path $peptide $weight"

# Optionally, attach to the tmux session (remove the comment below to enable)
# tmux attach-session -t "$session_name"
EOF

        # Make the batch script executable
        chmod +x "$batch_script"
        
        # Confirm the batch script creation
        echo "Tmux script created: $batch_script"

        # Submit the batch script to SLURM
        sbatch "$batch_script" && \
        echo "Submitted job: $job_name" || \
        echo "ERROR: Job submission failed for $job_name"
    done
done

echo "All tmux scripts created in the batch_scripts directory."

