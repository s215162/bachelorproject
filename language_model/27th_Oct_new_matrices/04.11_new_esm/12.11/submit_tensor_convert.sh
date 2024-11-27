#!/bin/bash
#SBATCH --job-name=process_tcr_vectors      # Job name
#SBATCH --output=process_tcr_vectors.out    # Standard output file
#SBATCH --error=process_tcr_vectors.err     # Standard error file
#SBATCH --time=01:00:00                     # Maximum run time (hh:mm:ss)
#SBATCH --mem=4G                            # Memory allocation (adjust as needed)
#SBATCH --cpus-per-task=1                  # Number of CPU cores


source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

# Run the Python script
python3 convert_tensor.py

