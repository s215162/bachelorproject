#!/bin/bash
#SBATCH --time=02:00:00             # Time limit (2 hours)
#SBATCH --mem=8G                   # Memory per job (24GB)
#SBATCH --array=0-47                # Array of jobs, each for different peptide and weight combination


# Activate Conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run Python script with peptide and weight
python3 esm_RAK_tcrb.py "RAKFKQLL" "TCRb"

# Deactivate environment after execution
conda deactivate
