#!/bin/bash

# Define arrays for peptides and partitions
peptides=("ELAGIGILTV" "GILGFVFTL" "GLCTLVAML" "LLWNGPMAV" "RAKFKQLL" "YLQPRTFLL")  
partitions=(0 1 2 3 4)
type=("with_swaps" "binders_only")


# Loop through each peptide
for peptide in "${peptides[@]}"; do
    # Get the peptide abbreviation (adjust as necessary)
    peptide_abbreviation=$(echo "$peptide" | cut -c 1-3)  # Adjust this line for your abbreviation logic
    
    # Loop through each partition
    for number in "${partitions[@]}"; do
        for type_used in "${type[@]}"; do
            # Call your scripts with the constructed paths and options
            ./weighed_full.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/${type_used}/partition${number}/partition_${number}_${peptide}" -o "AUC_${peptide_abbreviation}_${type_used}" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

            ./unweighed_full.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/${type_used}/partition${number}/partition_${number}_${peptide}" -o "AUC_${peptide_abbreviation}_${type_used}" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

            ./cdr3.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/${type_used}/partition${number}/partition_${number}_${peptide}" -o "AUC_${peptide_abbreviation}_${type_used}" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

            ./cdr3_alpha.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/${type_used}/partition${number}/partition_${number}_${peptide}" -o "AUC_${peptide_abbreviation}_${type_used}" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

            ./cdr3_beta.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/${type_used}/partition${number}/partition_${number}_${peptide}" -o "AUC_${peptide_abbreviation}_${type_used}" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition
        done
    done
done


