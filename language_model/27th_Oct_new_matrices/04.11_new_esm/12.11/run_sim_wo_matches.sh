#!/bin/bash
#SBATCH --job-name=get_similarity    # Job name
#SBATCH --output=output_%j.log       # Output log file (%j adds job ID)
#SBATCH --error=error_%j.log         # Error log file
#SBATCH --time=48:00:00              # Time limit (hh:mm:ss)
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=4            # CPUs per task
#SBATCH --mem=20G                     # Memory (adjust as needed)

# Activate the conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python3 sim_wo_dupes.py
