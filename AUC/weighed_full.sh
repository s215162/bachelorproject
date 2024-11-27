#!/bin/bash
# Here either /home/projects/vaccine/ or /home/ depending on computerome or HTC
start_time=$(date +%s)

# Paths
C2PATH="/home/projects/vaccine/people/morni/bin/tbcr_align"
HTCPATH="/home/people/morni/bin/tbcr_align"
CHAINS=("A1" "A2" "A3" "B1" "B2" "B3")  # Default chains
LABELCOL="peptide"
EXTRACOLS=("original_peptide" "binder" "partition" "original_index")
TBCRALIGN="$HTCPATH"
OUTPUTDIRECTORY=""

# Gets the argument from the command line (-f, -c, -s etc.)
while getopts ":f:c:s:l:e:o:" opt; do
  case ${opt} in
    f )
      INPUTFILE=$OPTARG # Saves the input filename
      ;;
    o )
      OUTPUTDIRECTORY="$OPTARG" # Get the output folder name
      echo "HERE ${OUTPUTDIRECTORY}"
      ;;
    c )
      # If -c is used, override the default chains
      CHAINS=("$OPTARG")  # Add the first option after -c
      while [[ ${!OPTIND} =~ ^[^-] ]]; do
        CHAINS+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    s )
      # Sets the path for tbcralign depending on which server you are using
      if [ "$OPTARG" == "htc" ]; then
        TBCRALIGN="$HTCPATH"
      elif [ "$OPTARG" == "c2" ]; then
        TBCRALIGN="$C2PATH"
      fi
      ;;
    e )
      # Extra cols are the column you want to be added to the final similarity output. 
      EXTRACOLS=("$OPTARG")  # Add the first option after -c
      while [[ ${!OPTIND} =~ ^[^-] ]]; do
        EXTRACOLS+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    l )
      LABELCOL="$OPTARG"  # Add the first option after -c
      ;;
    \? )
      echo "Usage: $0 -f <INPUTFILE> -o <OUTPUTDIRECTORY> -c <CHAINS> (ex: A1 A2 A3 B1 B2 B3) -s <SERVER> (c2/htc) -l <LABELCOL> -e <EXTRACOLS>"
      exit 1
      ;;
    : )
      echo "Invalid option: -$OPTARG requires an argument"
      exit 1
      ;;
  esac
done

# Shift the processed options so that $1, $2, etc. now refer to non-option arguments
shift $((OPTIND - 1))

# Creating the final output directory
OUTDIR="$(realpath "$(pwd)/output/${OUTPUTDIRECTORY}/")/"
mkdir -pv $OUTDIR

# Extract the basename without the extension using parameter expansion
filename=$(basename "$INPUTFILE")
basename_without_extension="${filename%.*}"
tmppath="${OUTDIR}${basename_without_extension}_TMP.txt" # Saves the intermediate TSV file for TBCRalign
tbcrtmp="${OUTDIR}${basename_without_extension}_TBCRraw_TMP.txt" # Saves the output of the TBCRalign
output_name="${OUTDIR}${basename_without_extension}_similarity_weighed.csv" # Final output name to save the similarity results

WEIGHTS="1,1,4,1,1,4"
# Variable for weight sum; change this value as needed
weight_sum=12

# Get the chains and extra columns as list to be used in the python part
string_chains=$(printf '%s' "$(IFS=','; echo "\"${CHAINS[*]}\"")")
string_extracols=$(printf '%s' "$(IFS=','; echo "\"${EXTRACOLS[*]}\"")")


# FILES MUST BE IN SAVED IN THE RIGHT FORMAT ; see below
# Embedded Python code to save the input df in the required format
/usr/bin/python3 -W ignore <<EOF
print('Saving tmp file')
import pandas as pd
filepath = "${INPUTFILE}"
tmppath = "${tmppath}"
chains = eval('${string_chains}').split(',')
df = pd.read_csv(filepath)
c='binder' if 'binder' in df.columns else 'target' if 'target' in df.columns else None
if c is None:
  c='binder'
  df[c] = 0.5 # Create a fake column
df[chains+[c]].to_csv(tmppath, sep='\t', header=False)
EOF

# Run TBCRalign
$TBCRALIGN -a -w $WEIGHTS $tmppath > $tbcrtmp

# Embedded Python code to parse and save results in the desired format
python3 -W ignore <<EOF
print('Reading TBCRalign output and saving selected columns')
import pandas as pd

# Read TBCRalign raw output
tbcr = pd.read_csv("${tbcrtmp}", comment='#', sep=' ', header=None, index_col=None)

# Define the relevant chains and extra columns
chains = eval('${string_chains}').split(',')
extra_cols = eval('${string_extracols}').split(',')
extra_cols = [x for x in extra_cols if len(x)>0]  # filter extra columns

# Set the column names for the TBCRalign output
qcols = [f'q_{x}' for x in chains + ['binder']]
dbcols = [f'db_{x}' for x in chains]
tbcr.columns = ['kind', 'source', 'q_index'] + qcols + ['target', 'db_index'] + dbcols + ['score'] + ['db_binder']

# We only want 'q_index', 'similarity', 'peptide', 'partition', and 'binder'
result = tbcr[['q_index', 'score', 'db_binder']].copy()

# Rename the columns as needed
result.columns = ['q_index', 'similarity', 'binder']

# Add peptide and partition columns from the original input file
original_df = pd.read_csv("${INPUTFILE}")
peptides = original_df['peptide']
partitions = original_df['partition']

# Match the peptide and partition values to the results
result['peptide'] = peptides
result['partition'] = partitions

# Divide similarity scores by weight_sum
result['similarity'] /= ${weight_sum}

# Save the result to a CSV file
output_filename = "${output_name}"
result.to_csv(output_filename, index=False)
print(f'Similarity data saved to {output_filename}')
EOF

# Cleaning up temporary files
rm ${OUTDIR}*TMP*.txt

# Record the end time
end_time=$(date +%s)

# Calculate the duration in seconds
duration=$(( end_time - start_time ))

# Convert the duration to HH:mm:ss format
hours=$(( duration / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$(( duration % 60 ))

# Format the output to HH:mm:ss
printf -v elapsed_time "%02d:%02d:%02d" $hours $minutes $seconds

# Output the elapsed time
echo "TBCRalign : Time taken: $elapsed_time"

