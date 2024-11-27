#!/bin/bash

start_time=$(date +%s)

HTCPATH="/home/people/morni/bin/tbcr_align"
CHAINS=("A1" "A2" "A3" "B1" "B2" "B3") # Default
LABELCOL="peptide"
EXTRACOLS=("original_peptide" "binder" "partition" "original_index")
TBCRALIGN="$HTCPATH"
OUTPUTDIRECTORY=""
SUM_WEIGHTS=0
WEIGHTS=""

# Gets the argument from the command line (-f, -c, -s etc.)
while getopts ":f:c:s:l:e:o:w:" opt; do
  case ${opt} in
    f )
      INPUTFILE=$OPTARG # Saves the inputfilename
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
      # Extra cols are the column you want to be added to the final distance matrix csv. For example it would be useful to keep the information of the "peptide", "binder", "original_index"
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
      echo "Usage: $0 -f <INPUTFILE> -o <OUTPUTDIRECTORY> -c <CHAINS> (ex: A1 A2 A3 B1 B2 B3) -s <SERVER> (c2/htc) -l <LABELCOL> -e <EXTRACOLS> -w <WEIGHTED>"
      exit 1
      ;;
    : )
      echo "Invalid option: -$OPTARG requires an argument"
      exit 1
      ;;
  esac
done


WEIGHTS="1,1,4,1,1,4"
SUM_WEIGHTS=12
# Shift the processed options so that $1, $2, etc. now refer to non-option arguments ; Forget this for now
shift $((OPTIND - 1))

echo "Weights: $WEIGHTS"
echo "Sum of Weights: $SUM_WEIGHTS"

# Get the chains and extracols as list to be used in the python part
string_chains=$(printf '%s' "$(IFS=','; echo "\"${CHAINS[*]}\"")")
string_extracols=$(printf '%s' "$(IFS=','; echo "\"${EXTRACOLS[*]}\"")")

# Remove the quotes and replace commas with underscores
output_filename_chains=$(echo "$string_chains" | tr -d '"' | tr ',' '_')

OUTDIR="$(realpath "$(pwd)/output/${OUTPUTDIRECTORY}/")/"
mkdir -pv $OUTDIR

# Extract the basename without the extension using parameter expansion
filename=$(basename "$INPUTFILE")
basename_without_extension="${filename%.*}"
tmppath="${OUTDIR}${basename_without_extension}_TMP.txt" # Saves the intermediate tsv file for TBCRalign
tbcrtmp="${OUTDIR}${basename_without_extension}_TBCRraw_TMP.txt" # Saves the output of the TBCRalign to be read to convert to distmatrix
output_name="${OUTDIR}${basename_without_extension}_${output_filename_chains}_TBCR_distmatrix.csv" # final output name to save the dm

# FILES MUST BE IN SAVED IN THE RIGHT FORMAT ; see below
  # Embedded python code to save the input df in the required format
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


# Embedded Python code to recover the distance matrix from the AvA distance list
python3 -W ignore <<EOF
print('Reading TBCRalign output and converting to dist_matrix')
import pandas as pd
import numpy as np
# Accessing the Bash variable in Python
output_filename = "${output_name}"
print('Read TBCR raw')
tbcr = pd.read_csv("${tbcrtmp}", comment='#', sep =' ', header=None, index_col=None)

tbcr.to_csv(output_filename)
print(f'TBCR distmatrix saved at {output_filename}')
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

