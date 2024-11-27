#!/bin/bash
#SBATCH --job-name=LLWNGPMAV_TCRa_swaps
#SBATCH --output="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/batch_scripts/LLWNGPMAV_TCRa_swaps_%j.out"
#SBATCH --error="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/batch_scripts/LLWNGPMAV_TCRa_swaps_%j.err"
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G

# Activate conda environment
source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python3 /net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/run_esm_swaps.py LLWNGPMAV TCRa

# Deactivate the conda environment
conda deactivate
