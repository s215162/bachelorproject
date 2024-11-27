#!/bin/bash
#SBATCH --job-name="LLWNGPMAV_TCRa_binders"
#SBATCH --output="LLWNGPMAV_TCRa_binders_%j.out"
#SBATCH --error="LLWNGPMAV_TCRa_binders_%j.err"
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --constraint=gpu
#SBATCH --mem=16G

# Activate conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python /net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/run_esm_binders_new.py LLWNGPMAV TCRa

# Deactivate the conda environment
conda deactivate
