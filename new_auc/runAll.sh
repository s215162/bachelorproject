#! /bin/bash

echo "ELA starting"
./plswork.sh -f "/home/projects2/emison/auc/divided_data/Divided_Data/with_swaps/partition0/partition_0_ELAGIGILTV" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "GIL starting"
./plswork.sh -f "/home/projects2/emison/new_auc/training_data_by_peptide/training_GILGFVFTL_data.csv" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "GLC starting"
./plswork.sh -f "/home/projects2/emison/new_auc/training_data_by_peptide/training_GLCTLVAML_data.csv" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "LLW starting"

./plswork.sh -f "/home/projects2/emison/new_auc/training_data_by_peptide/training_LLWNGPMAV_data.csv" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "RAK starting"

./plswork.sh -f "/home/projects2/emison/new_auc/training_data_by_peptide/training_RAKFKQLL_data.csv" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "YLQ starting"

./plswork.sh -f "/home/projects2/emison/new_auc/training_data_by_peptide/training_YLQPRTFLL_data.csv" -o "sim_matrix" -c A1 A2 A3 B1 B2 B3 -s htc -l binder -e peptide binder partition

echo "finished"

