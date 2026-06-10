#!/usr/bin/bash

# Directory to read files from (can be passed as an argument, or use current directory)
REPCALC_DIR="/mnt/md0/Projects/MonsonLab/MAV/vdjserver/d2296f11-4358-4972-9abd-263dde5c3287-007"
IGBLAST_DIR="/mnt/md0/Projects/MonsonLab/MAV/vdjserver/50524118-e81a-47f6-bc9d-2345969b9bfe-007"
SUMMARY_DIR="/mnt/md0/Projects/MonsonLab/MAV/summary"

# Check if directory exists
if [ ! -d "$REPCALC_DIR" ]; then
    echo "Directory does not exist: $REPCALC_DIR"
    exit 1
fi

# Collect prefixes
prefixes=()

while IFS= read -r file; do
    filename=$(basename "$file")
    prefix="${filename%%.*}"
    prefixes+=("$prefix")
done < <(find "$REPCALC_DIR" -maxdepth 1 -type f)

# Output unique prefixes
unique_prefixes=($(printf "%s\n" "${prefixes[@]}" | sort -u))

for prefix in "${unique_prefixes[@]}"; do
    echo "Processing repertoire: ${prefix}"

    MUTATION_FILE="${REPCALC_DIR}/${prefix}.igblast.makedb.gene.clone.mutations.airr.tsv"
    if [ -f "${MUTATION_FILE}.gz" ]; then
        echo "Found gzipped file: ${MUTATION_FILE}.gz, unzipping..."
        gunzip -f "${MUTATION_FILE}.gz"
    fi

    IGBLAST_FILE="${IGBLAST_DIR}/${prefix}.igblast.airr.tsv"
    if [ -f "${IGBLAST_FILE}.gz" ]; then
        echo "Found gzipped file: ${IGBLAST_FILE}.gz, unzipping..."
        gunzip -f "${IGBLAST_FILE}.gz"
    fi

    AA_PROP_FILE="${REPCALC_DIR}/${prefix}.igblast.makedb.gene.clone.aa_properties.airr.tsv"
    if [ -f "${AA_PROP_FILE}.gz" ]; then
        echo "Found gzipped file: ${AA_PROP_FILE}.gz, unzipping..."
        gunzip -f "${AA_PROP_FILE}.gz"
    fi
    
    python3 /mnt/md0/s236922/cowell-lab/monson-lab-projects/RHAB/rhab_summary_v5.py ${MUTATION_FILE} ${IGBLAST_FILE} ${AA_PROP_FILE} ${SUMMARY_DIR}/${prefix}.summary.tsv

done

# set -eu pipefail

# REPCALC_DIR="/mnt/md0/Projects/MonsonLab/MAV/vdjserver/d2296f11-4358-4972-9abd-263dde5c3287-007"
# IGBLAST_DIR="/mnt/md0/Projects/MonsonLab/MAV/vdjserver/50524118-e81a-47f6-bc9d-2345969b9bfe-007"
# SUMMARY_DIR="/mnt/md0/Projects/MonsonLab/MAV/summary"
# REPERTOIRE_JSON="/mnt/md0/Projects/MonsonLab/MAV/repertoires.airr.json"

# # Extract prefix (repertoire_id) and sequencing filename
# mapfile -t entries < <(
#     jq -r '
#         .Repertoire[] |
#         .repertoire_id as $rid |
#         .sample[0].sequencing_files.filename as $fname |
#         "\($rid)\t\($fname)"
#     ' "$REPERTOIRE_JSON"
# )

# for entry in "${entries[@]}"; do
#     prefix=$(echo "$entry" | cut -f1)
#     fastq_filename=$(echo "$entry" | cut -f2)

#     echo "Processing: $prefix"
#     echo " FASTQ filename: $fastq_filename"

#     # Convert FASTQ → IGBlast TSV filename
#     igblast_file="${fastq_filename%.fastq.gz}.igblast.airr.tsv"

#     python3 /mnt/md0/s236922/cowell-lab/monson-lab-projects/RHAB/rhab_summary_v5.py \
#         "${REPCALC_DIR}/${prefix}.igblast.makedb.gene.clone.mutations.airr.tsv" \
#         "${IGBLAST_DIR}/${igblast_file}" \
#         "${REPCALC_DIR}/${prefix}.igblast.makedb.gene.clone.aa_properties.airr.tsv" \
#         "${SUMMARY_DIR}/${prefix}.summary.tsv"
# done
