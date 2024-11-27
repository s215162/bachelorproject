#!/bin/bash
#SBATCH --time=02:00:00             # Time limit (2 hours)
#SBATCH --mem=8G                   # Memory per job (24GB)
#SBATCH --array=0-47                # Array of jobs, each for different peptide and weight combination


# Activate Conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run Python script with peptide and weight
python3 get_binders.py "ELAGIGILTV" 
python3 get_binders.py "GILGFVFTL" 
python3 get_binders.py "GLCTLVAML" 
python3 get_binders.py "LLWNGPMAV" 
python3 get_binders.py "RAKFKQLL" 
python3 get_binders.py "YLQPRTFLL" 


# Deactivate environment after execution
conda deactivate

echo "Job for Peptide: $peptide" >> logs/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out

