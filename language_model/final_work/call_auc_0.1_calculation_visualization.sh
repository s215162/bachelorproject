#!/bin/bash
#SBATCH --job-name=auc_0.1
#SBATCH --output=logs/%x_%j.out    # Log file (one per job)
#SBATCH --error=logs/%x_%j.err     # Error log
#SBATCH --time=1-00:00:00          # Max runtime (1 days)
#SBATCH --mem=128G                 # Memory per job

# Activate environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script for the given peptide
python3 /net/mimer/mnt/tank/projects2/emison/language_model/final_work/auc_0.1_calculation_visualization.py
