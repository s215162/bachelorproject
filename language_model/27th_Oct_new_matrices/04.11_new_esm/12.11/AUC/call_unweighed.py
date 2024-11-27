#!/bin/bash
#SBATCH --time=02:00:00             # Time limit (2 hours)
#SBATCH --mem=8G                   # Memory per job (24GB)


# Activate Conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run Python script with peptide and weight
python3 unweighed.py "ELAGIGILTV"
python3 unweighed.py "GILGFVFTL"
python3 unweighed.py "GLCTLVAML"
python3 unweighed.py "LLWNGPMAV"
python3 unweighed.py "RAKFKQLL"
python3 unweighed.py "YLQPRTFLL"


# Deactivate environment after execution
conda deactivate
