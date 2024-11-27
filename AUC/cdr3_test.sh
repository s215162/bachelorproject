
#!/bin/bash
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
      INPUTFILE=$OPTARG
      ;;
    o )
      OUTPUTDIRECTORY="$OPTARG"
      ;;
    c )
      CHAINS=("$OPTARG")
      while [[ ${!OPTIND} =~ ^[^-] ]]; do
        CHAINS+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    s )
      if [ "$OPTARG" == "htc" ]; then
        TBCRALIGN="$HTCPATH"
      elif [ "$OPTARG" == "c2" ]; then
        TBCRALIGN="$C2PATH"
      fi
      ;;
    e )
      EXTRACOLS=("$OPTARG")
      while [[ ${!OPTIND} =~ ^[^-] ]]; do
        EXTRACOLS+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    l )
      LABELCOL="$OPTARG"
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

shift $((OPTIND - 1))

OUTDIR="$(realpath "$(pwd)/output/${OUTPUTDIRECTORY}/")/"
mkdir -pv $OUTDIR

filename=$(basename "$INPUTFILE")
basename_without_extension="${filename%.*}"
tmppath="${OUTDIR}${basename_without_extension}_TMP.txt"
tbcrtmp="${OUTDIR}${basename_without_extension}_TBCRraw_TMP.txt"
output_name="${OUTDIR}${basename_without_extension}_similarity_CDR3_ALPHA.csv"

WEIGHTS="0,0,1,0,0,0"
weight_sum=1

string_chains=$(printf '%s' "$(IFS=','; echo "\"${CHAINS[*]}\"")")
string_extracols=$(printf '%s' "$(IFS=','; echo "\"${EXTRACOLS[*]}\"")")

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

$TBCRALIGN -a -w $WEIGHTS $tmppath > $tbcrtmp

python3 -W ignore <<EOF
print('Reading TBCRalign output and saving selected columns')
import pandas as pd

tbcr = pd.read_csv("${tbcrtmp}", comment='#', sep=' ', header=None, index_col=None)

chains = eval('${string_chains}').split(',')
extra_cols = eval('${string_extracols}').split(',')
extra_cols = [x for x in extra_cols if len(x)>0]

qcols = [f'q_{x}' for x in chains + ['binder']]
dbcols = [f'db_{x}' for x in chains]
tbcr.columns = ['kind', 'source', 'q_index'] + qcols + ['target', 'db_index'] + dbcols + ['score'] + ['db_binder']

result = tbcr[['score', 'db_binder']].copy()

result.columns = ['similarity', 'binder']

result['similarity'] /= ${weight_sum}

output_filename = "${output_name}"
result.to_csv(output_filename, index=False)
print(f'Similarity data saved to {output_filename}')
EOF

rm ${OUTDIR}*TMP*.txt

end_time=$(date +%s)
duration=$(( end_time - start_time ))
hours=$(( duration / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$(( duration % 60 ))

printf -v elapsed_time "%02d:%02d:%02d" $hours $minutes $seconds
echo "TBCRalign : Time taken: $elapsed_time"

