#!/bin/bash
#SBATCH --job-name=emb_to_sim  # Job name for easy identification in SLURM queue
#SBATCH --output=/home/projects2/emison/logs/emb_to_sim_%j.log  # Path for stdout logs; %j is job ID
#SBATCH --error=/home/projects2/emison/logs/emb_to_sim_%j.err   # Path for stderr logs; %j is job ID
#SBATCH --time=00:20:00                                         # Max wall time limit; adjust as necessary
#SBATCH --partition=cpu                                     # SLURM partition to use; adjust as needed
#SBATCH --nodes=1                                               # Number of nodes; 1 node for this job
#SBATCH --ntasks=1                                              # Number of tasks (processes) per node
#SBATCH --cpus-per-task=1                                       # Number of CPUs per task; adjust as needed
#SBATCH --mem=8G                                                # Memory allocation per task; adjust as needed

# Load Conda
source /home/projects2/emison/language_model/miniconda3/bin/activate

# Activate the specific environment
conda activate esm2_emison

# Set output directory for saving the generated similarity matrices
OUTPUT_DIR="/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/cos_sims/"

# Array of commands to run
commands=(
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type alpha --sequence cdr1a --output_dis cdr1a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type alpha --sequence cdr2a --output_dis cdr2a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type alpha --sequence cdr3a --output_dis cdr3a --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type alpha --sequence tcra --output_dis tcra --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide ELAGIGILTV --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GILGFVFTL --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide GLCTLVAML --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide LLWNGPMAV --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide RAKFKQLL --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type beta --sequence cdr1b --output_dis cdr1b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type beta --sequence cdr2b --output_dis cdr2b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type beta --sequence cdr3b --output_dis cdr3b --output_dir $OUTPUT_DIR"
  "python3 emb_to_sim.py --peptide YLQPRTFLL --type beta --sequence tcrb --output_dis tcrb --output_dir $OUTPUT_DIR"
)

# Loop over each command in the commands array
# Each command will be submitted as a separate SLURM job
for cmd in "${commands[@]}"
do
  echo "Submitting job for: $cmd"    # Print the command being executed (for logging purposes)
  echo -e "#!/bin/bash\n$cmd" | sbatch  # Submit the command as an independent SLURM job
done

# Deactivate Conda environment
conda deactivate

echo "All tasks submitted."  # Message indicating all commands have been submitted

