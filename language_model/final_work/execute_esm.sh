#!/bin/bash
peptides=(
    "ELAGIGILTV"
    "GILGFVFTL"
    "GLCTLVAML"
    "LLWNGPMAV"
    "RAKFKQLL"
    "YLQPRTFLL"
)
for peptide in "${peptides[@]}"; do
    sbatch --wrap="python /net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_matrix_full_sequence.py $peptide"
    sbatch --wrap="python /net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_mean_full_sequence.py $peptide"
    sbatch --wrap="python /net/mimer/mnt/tank/projects2/emison/language_model/final_work/esm_sum_full_sequence.py $peptide"
done

