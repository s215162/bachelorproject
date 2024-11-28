#!/bin/bash
#SBATCH --time=02:00:00             # Time limit (2 hours)
#SBATCH --mem=16G                   # Memory per job (8GB)


# Activate Conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run Python script with peptide and weight
python3 DBSCAN.py

# Deactivate environment after execution
conda deactivate