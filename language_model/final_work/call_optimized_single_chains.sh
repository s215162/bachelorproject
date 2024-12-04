#!/bin/bash

# List of peptides to process
peptides=(
    "ELAGIGILTV"
    "GILGFVFTL"
    "GLCTLVAML"
    "LLWNGPMAV"
    "RAKFKQLL"
    "YLQPRTFLL"
)

# Path to the Python script
sequence_script="/net/mimer/mnt/tank/projects2/emison/language_model/final_work/optimized_single_chains.py"

# Loop through peptides and submit jobs
for peptide in "${peptides[@]}"; do
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=esm_${peptide}
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=128G

source /home/projects2/emison/language_model/miniconda3/bin/activate
conda activate esm2_emison

python3 $sequence_script "$peptide" "TCRb"
EOF
done


# change TCRb to TCRa depending on whether you want to embed the alpha or beta chain