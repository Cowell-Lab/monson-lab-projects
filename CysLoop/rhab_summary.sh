 #!/bin/bash

# Directory to read files from (can be passed as an argument, or use current directory)
DIR="/Users/s236922/code/cowell-lab/monson-lab-projects/CysLoop/jobs/repcalc/0d895ccd-9649-4a56-8a13-c1591a148fca-007/0d895ccd-9649-4a56-8a13-c1591a148fca-007"
IGBLAST_DIR="/Users/s236922/code/cowell-lab/monson-lab-projects/CysLoop/jobs/igblast/3ec8bcc8-db20-4a8e-965e-e8587914dca8-007/3ec8bcc8-db20-4a8e-965e-e8587914dca8-007"
SOLO_BLAST_DIR="/Users/s236922/code/cowell-lab/monson-lab-projects/CysLoop/jobs/igblast/942c6f01-606a-4307-bf41-e63eb2cd77c4-007/942c6f01-606a-4307-bf41-e63eb2cd77c4-007"
SUMMARY_DIR="/Users/s236922/code/cowell-lab/monson-lab-projects/CysLoop/summary"

# Check if directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory does not exist: $DIR"
    exit 1
fi

# Collect prefixes
prefixes=()

while IFS= read -r file; do
    filename=$(basename "$file")
    prefix="${filename%%.*}"
    prefixes+=("$prefix")
done < <(find "$DIR" -maxdepth 1 -type f)

# Output unique prefixes
unique_prefixes=($(printf "%s\n" "${prefixes[@]}" | sort -u))

for prefix in "${unique_prefixes[@]}"; do
    echo "Processing repertoire: ${prefix}"
    
    # gunzip ${DIR}/${prefix}.igblast.makedb.gene.clone.mutations.airr.tsv.gz
    # gunzip ${DIR}/${prefix}.igblast.makedb.gene.clone.aa_properties.airr.tsv.gz

    if [[ "$prefix" != "1bf67ffb-75de-40b3-b33c-580f013f1b41" ]]; then
        # gunzip ${IGBLAST_DIR}/${prefix}.igblast.airr.tsv.gz
        python3 ~/code/cowell-lab/monson-lab-projects/RHAB/rhab_summary_v5.py ${DIR}/${prefix}.igblast.makedb.gene.clone.mutations.airr.tsv ${IGBLAST_DIR}/${prefix}.igblast.airr.tsv ${DIR}/${prefix}.igblast.makedb.gene.clone.aa_properties.airr.tsv ${SUMMARY_DIR}/${prefix}.summary.tsv
    else
        # gunzip ${SOLO_BLAST_DIR}/${prefix}.igblast.airr.tsv.gz
        python3 ~/code/cowell-lab/monson-lab-projects/RHAB/rhab_summary_v5.py ${DIR}/${prefix}.igblast.makedb.gene.clone.mutations.airr.tsv ${SOLO_BLAST_DIR}/${prefix}.igblast.airr.tsv ${DIR}/${prefix}.igblast.makedb.gene.clone.aa_properties.airr.tsv ${SUMMARY_DIR}/${prefix}.summary.tsv
    fi
done
