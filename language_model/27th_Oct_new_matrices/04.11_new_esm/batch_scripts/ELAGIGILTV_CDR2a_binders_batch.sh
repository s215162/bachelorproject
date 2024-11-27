#!/bin/bash
#SBATCH --job-name="ELAGIGILTV_CDR2a_binders"
#SBATCH --output="ELAGIGILTV_CDR2a_binders_%j.out"
#SBATCH --error="ELAGIGILTV_CDR2a_binders_%j.err"
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --constraint=gpu
#SBATCH --mem=16G

# Activate conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python /net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/run_esm_binders_new.py ELAGIGILTV CDR2a

# Deactivate the conda environment
conda deactivate
