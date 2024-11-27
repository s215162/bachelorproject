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

# Set the directory to save the batch scripts
batch_dir="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/batch_scripts"

# Ensure the batch script directory exists
mkdir -p "$batch_dir"

# Activate the conda environment for esm2
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Script path for the Python script that will be executed
script_path="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/run_esm_binders_new.py"

# Iterate over each peptide and weight to submit a job to SLURM
for peptide in "${peptides[@]}"; do
    for weight in "${weights[@]}"; do
        # Create a unique job name for each peptide-weight combination
        job_name="${peptide}_${weight}_binders"

        # Create a temporary SLURM batch script for this job in the specified directory
        batch_script="${batch_dir}/${job_name}_batch.sh"
        
        # Write the SLURM script to the batch script file
        cat <<EOF > "$batch_script"
#!/bin/bash
#SBATCH --job-name="$job_name"
#SBATCH --output="${job_name}_%j.out"
#SBATCH --error="${job_name}_%j.err"
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --constraint=gpu
#SBATCH --mem=16G

# Activate conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python $script_path $peptide $weight

# Deactivate the conda environment
conda deactivate
EOF

        # Confirm the batch script creation
        echo "Batch script created: $batch_script"

        # Submit the batch script to SLURM
        sbatch "$batch_script" && \
        echo "Submitted job: $job_name" || \
        echo "ERROR: Job submission failed for $job_name"

    done
done

# Wait for all jobs to finish
wait

# Final confirmation message after all jobs have been submitted
echo "All jobs submitted."

