#!/bin/bash

# Create a new tmux session named 'embeddings_sim'
tmux new-session -d -s embeddings_sim

# Function to run commands in a new tmux window
run_in_tmux() {
    tmux new-window -t embeddings_sim -n "$1"
    tmux send-keys -t embeddings_sim:"$1" "python3 emb_to_sim.py --peptide \"$2\" --type \"$3\" --sequence \"$4\" --output_dis \"$5\" --output_dir \"$6\"" C-m
}

# Run commands in separate tmux windows
run_in_tmux "beta_full_GILGFVFTL" "GILGFVFTL" "beta" "full" "full_GILGFVFTL" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/full/GILGFVFTL/"
run_in_tmux "beta_cdr1_GILGFVFTL" "GILGFVFTL" "beta" "cdr1" "cdr1_GILGFVFTL" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/beta/cdr1/GILGFVFTL/"
run_in_tmux "alpha_cdr1_LLWNGPMAV" "LLWNGPMAV" "alpha" "cdr1" "cdr1_LLWNGPMAV" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr1/LLWNGPMAV/"
run_in_tmux "alpha_cdr2_LLWNGPMAV" "LLWNGPMAV" "alpha" "cdr2" "cdr2_LLWNGPMAV" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr2/LLWNGPMAV/"
run_in_tmux "alpha_cdr3_LLWNGPMAV" "LLWNGPMAV" "alpha" "cdr3" "cdr3_LLWNGPMAV" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/cdr3/LLWNGPMAV/"
run_in_tmux "alpha_full_LLWNGPMAV" "LLWNGPMAV" "alpha" "full" "full_LLWNGPMAV" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/full/LLWNGPMAV/"
run_in_tmux "alpha_full_GLCTLVAML" "GLCTLVAML" "alpha" "full" "full_GLCTLVAML" "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/alpha/full/GLCTLVAML/"

# Attach to the tmux session
tmux attach-session -t embeddings_sim

