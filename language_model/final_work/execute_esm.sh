#!/bin/bash
#SBATCH --job-name=new_esm_matrix
#SBATCH --output=logs/%x_%j.out    # Log file (one per job)
#SBATCH --error=logs/%x_%j.err     # Error log
#SBATCH --time=3-00:00:00          # Max runtime (7 days)
#SBATCH --mem=128G                 # Memory per job

# List of peptides to process
peptides=(
    "ELAGIGILTV"
    "GILGFVFTL"
    "GLCTLVAML"
    "LLWNGPMAV"
    "RAKFKQLL"
    "YLQPRTFLL"
)

# Paths to the Python scripts
matrix_script="/net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_matrix_full_sequence.py"
mean_script="/net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_mean_full_sequence.py"
sum_script="/net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_sum_full_sequence.py"

# Loop through peptides and submit jobs for each one
for peptide in "${peptides[@]}"; do
    # Submit a job for each peptide
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=esm_matrix_${peptide}
#SBATCH --output=logs/${peptide}_%j.out
#SBATCH --error=logs/${peptide}_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=128G

# Activate conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the first Python script for the given peptide
#python3 $mean_script "$peptide"

# Run the second Python script for the given peptide
#python3 $sum_script "$peptide"

# Run the third Python script for the given peptide
python3 $matrix_script "$peptide"
EOF
done