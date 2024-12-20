#!/bin/bash

start_time=$(date +%s)

# Set paths here
HTCPATH=""  # Path to tbcr_align program
OUTPUTDIRECTORY=""  # Where to save processed data

# Default configurations
CHAINS=("A1" "A2" "A3" "B1" "B2" "B3")
LABELCOL="peptide"
EXTRACOLS=("original_peptide" "binder" "partition" "original_index")
TBCRALIGN="$HTCPATH"
WEIGHTED=""
SUM_WEIGHTS=0

# Process command line args
while getopts ":f:c:s:l:e:o:w:" opt; do
 case ${opt} in
   f ) INPUTFILE=$OPTARG ;;
   o ) OUTPUTDIRECTORY="$OPTARG"
       echo "Output directory: ${OUTPUTDIRECTORY}" ;;
   c ) # Override default chains
       CHAINS=("$OPTARG")
       while [[ ${!OPTIND} =~ ^[^-] ]]; do
         CHAINS+=("${!OPTIND}")
         OPTIND=$((OPTIND + 1))
       done ;;
   s ) # Set tbcralign path based on server
       if [ "$OPTARG" == "htc" ]; then
         TBCRALIGN="$HTCPATH"
       elif [ "$OPTARG" == "c2" ]; then
         TBCRALIGN="$C2PATH"
       fi ;;
   e ) # Add extra columns to final matrix
       EXTRACOLS=("$OPTARG")
       while [[ ${!OPTIND} =~ ^[^-] ]]; do
         EXTRACOLS+=("${!OPTIND}")
         OPTIND=$((OPTIND + 1))
       done ;;
   l ) LABELCOL="$OPTARG" ;;
   w ) WEIGHTED="$OPTARG" ;;
   \? )
       echo "Usage: $0 -f <INPUTFILE> -o <OUTPUTDIRECTORY> -c <CHAINS> -s <SERVER> -l <LABELCOL> -e <EXTRACOLS> -w <WEIGHTED>"
       exit 1 ;;
   : )
       echo "Invalid option: -$OPTARG requires an argument"
       exit 1 ;;
 esac
done

read -p "Is this weighted? (y/n): " weighted_response

# Set weights based on chain combinations
if [[ "${CHAINS[@]}" =~ "A1" && "${CHAINS[@]}" =~ "A2" && "${CHAINS[@]}" =~ "A3" && "${CHAINS[@]}" =~ "B1" && "${CHAINS[@]}" =~ "B2" && "${CHAINS[@]}" =~ "B3" ]]; then
   if [ "$weighted_response" == "y" ]; then
       WEIGHTS="1,1,4,1,1,4"
       SUM_WEIGHTS=12
   else
       WEIGHTS="1,1,1,1,1,1"
       SUM_WEIGHTS=6
   fi
elif [[ "${CHAINS[@]}" =~ "A1" && "${CHAINS[@]}" =~ "B1" ]]; then
   WEIGHTS="1,0,0,1,0,0"
   SUM_WEIGHTS=2
elif [[ "${CHAINS[@]}" =~ "A2" && "${CHAINS[@]}" =~ "B2" ]]; then
   WEIGHTS="0,1,0,0,1,0"
   SUM_WEIGHTS=2
elif [[ "${CHAINS[@]}" =~ "A3" && "${CHAINS[@]}" =~ "B3" ]]; then
   WEIGHTS="0,0,1,0,0,1"
   SUM_WEIGHTS=2
elif [[ "${CHAINS[@]}" =~ "A3" ]]; then
   WEIGHTS="0,0,1,0,0,0"
   SUM_WEIGHTS=1
elif [[ "${CHAINS[@]}" =~ "B3" ]]; then
   WEIGHTS="0,0,0,0,0,1"
   SUM_WEIGHTS=1
else
   echo "Unsupported chain combination. Exiting."
   exit 1
fi

WEIGHTS="1,1,4,1,1,4"

echo "Weights: $WEIGHTS"
echo "Sum of Weights: $SUM_WEIGHTS"

shift $((OPTIND - 1))

# Setup output paths
OUTDIR="$(realpath "$(pwd)/output/${OUTPUTDIRECTORY}/")/"
mkdir -pv $OUTDIR

# Prepare lists
string_chains=$(printf '%s' "$(IFS=','; echo "\"${CHAINS[*]}\"")")
string_extracols=$(printf '%s' "$(IFS=','; echo "\"${EXTRACOLS[*]}\"")")
output_filename_chains=$(echo "$string_chains" | tr -d '"' | tr ',' '_')

# filenames initiate
filename=$(basename "$INPUTFILE")
basename_without_extension="${filename%.*}"
tmppath="${OUTDIR}${basename_without_extension}_TMP.txt"
tbcrtmp="${OUTDIR}${basename_without_extension}_TBCRraw_TMP.txt"
output_name="${OUTDIR}${basename_without_extension}_binders_TBCR_distmatrix_FULL_weighed.csv"

# formatting and saving input data
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
 df[c] = 0.5
df[chains+[c]].to_csv(tmppath, sep='\t', header=False)
EOF

# Run TBCRalign
$TBCRALIGN -a -w $WEIGHTS $tmppath > $tbcrtmp

# Process results and create distance matrix
python3 -W ignore <<EOF
print('Converting TBCRalign output to distance matrix')
import pandas as pd
import numpy as np

output_filename = "${output_name}"
tbcr = pd.read_csv("${tbcrtmp}", comment='#', sep =' ', header=None, index_col=None)
chains = eval('${string_chains}').split(',')
extra_cols = eval('${string_extracols}').split(',')
extra_cols = [x for x in extra_cols if len(x)>0]

qcols = [f'q_{x}' for x in chains + ['binder']]
dbcols = [f'db_{x}' for x in chains]

print('Processing data...')
tbcr.columns = ['kind', 'source', 'q_index'] + qcols + ['target', 'db_index'] + dbcols + ['score'] + ['db_binder']

fake_tbcr = tbcr.copy()
fake_tbcr.columns = ['kind', 'name', 'db_index']+ dbcols + ['db_binder', 'target','q_index'] + [x for x in qcols if not x.endswith('binder')] + ['score','q_binder']

tbcr_cat = pd.concat([tbcr, fake_tbcr])
tbcr_idx_scores = tbcr_cat.drop(columns=[x for x in tbcr_cat.columns if 'index' not in x and 'score' not in x])

index_tcrs = tbcr_cat.drop_duplicates(subset=['q_index'])
index_tcrs['tcr'] = index_tcrs[[x for x in qcols if x != 'q_binder' and x != 'q_label' and x !='q_target']].sum(axis=1)
index_tcrs = index_tcrs[['q_index', 'tcr']]

tbcr_idx_scores = pd.concat([tbcr_idx_scores, pd.DataFrame([[x, x, ${SUM_WEIGHTS}.] for x in tbcr_idx_scores.q_index.unique()], columns=['q_index', 'db_index', 'score'])])
dist_matrix = 1-(tbcr_idx_scores.pivot_table(index='q_index', columns='db_index', values='score')/${SUM_WEIGHTS})

original_filename = "${INPUTFILE}"
original_df = pd.read_csv(original_filename)
original_df['tcr'] = original_df[chains].sum(axis=1)

len_before = len(original_df)
original_df = original_df.merge(index_tcrs[['q_index','tcr']], left_on='tcr', right_on='tcr').drop_duplicates(['q_index','tcr'])
len_after = len(original_df)
extra_cols = [x for x in extra_cols if x in original_df.columns]

print(f'Saving distance matrix with extra columns {extra_cols}')
if len_before!=len_after:
   print(f'Warning: length mismatch - before: {len_before}, after: {len_after}')

dist_matrix = dist_matrix.merge(original_df.set_index('q_index')[extra_cols], left_index=True, right_index=True)
dist_matrix.to_csv(output_filename)
print(f'Distance matrix saved to {output_filename}')
EOF

# Cleanup
rm ${OUTDIR}*TMP*.txt

# Calculate runtime
end_time=$(date +%s)
duration=$(( end_time - start_time ))
hours=$(( duration / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$(( duration % 60 ))
printf -v elapsed_time "%02d:%02d:%02d" $hours $minutes $seconds
echo "Runtime: $elapsed_time"