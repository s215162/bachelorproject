#!/bin/bash
#SBATCH --time=02:00:00             # Time limit (2 hours)
#SBATCH --mem=24G                   # Memory per job (24GB)
#SBATCH --array=0-47                # Array of jobs, each for different peptide and weight combination
#SBATCH --output=logs/job_%A_%a.out # Standard output (logs)
#SBATCH --error=logs/job_%A_%a.err  # Standard error (logs)

# Make sure logs directory exists
mkdir -p logs

# List of peptides and weights
peptides=("ELAGIGILTV" "GILGFVFTL" "GLCTLVAML" "LLWNGPMAV" "RAKFKQLL" "YLQPRTFLL")
weights=("CDR1a" "CDR1b" "CDR2a" "CDR2b" "CDR3a" "CDR3b" "TCRa" "TCRb")

# Compute peptide and weight indices based on task ID
peptide_index=$(( SLURM_ARRAY_TASK_ID / ${#weights[@]} ))
weight_index=$(( SLURM_ARRAY_TASK_ID % ${#weights[@]} ))

peptide=${peptides[$peptide_index]}
weight=${weights[$weight_index]}

# Print debug information to logs
echo "Task ID: $SLURM_ARRAY_TASK_ID" >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
echo "Peptide: $peptide" >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
echo "Weight: $weight" >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

# Activate Conda environment
echo "Activating Conda environment..." >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run Python script with peptide and weight
echo "Running Python script..." >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
python3 new_similarity.py "$peptide" "$weight" >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out 2>&1

# Deactivate environment after execution
echo "Deactivating Conda environment..." >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out
conda deactivate

echo "Job for Peptide: $peptide and Weight: $weight completed." >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

