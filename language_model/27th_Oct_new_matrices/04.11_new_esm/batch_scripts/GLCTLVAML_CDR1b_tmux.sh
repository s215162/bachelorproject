#!/bin/bash
# This script runs the Python script in a tmux session

# Start a new tmux session and run the Python script
tmux new-session -d -s "GLCTLVAML_CDR1b" "python /net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/new_esm/run_esm_new.py GLCTLVAML CDR1b"

# Optionally, attach to the tmux session (remove the comment below to enable)
# tmux attach-session -t "GLCTLVAML_CDR1b"
