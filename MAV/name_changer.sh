#!/usr/bin/bash
SUMMARY_DIR="/mnt/md0/Projects/MonsonLab/MAV/summary"
REPERTOIRE_JSON="/mnt/md0/Projects/MonsonLab/MAV/repertoires.airr.json"

mapfile -t ENTRIES < <(
    jq -r '
        .Repertoire[] |
        .repertoire_id as $RID |
        .sample[0].sequencing_files.filename as $FNAME |
        "\($RID)\t\($FNAME)"
    ' "$REPERTOIRE_JSON"
)

for ENTRY in "${ENTRIES[@]}"; do
    REPERTOIRE_ID=$(echo "$ENTRY" | cut -f1)
    SEQ_FILENAME=$(echo "$ENTRY" | cut -f2)
    SAMPLE_ID=$(echo "$SEQ_FILENAME" | cut -d . -f1)
    SUMMARY_FILENAME=$REPERTOIRE_ID".summary.tsv"

    echo "Processing: $REPERTOIRE_ID |  Summary filesname: $SUMMARY_FILENAME | Sample ID: $SAMPLE_ID"

    mv ${SUMMARY_DIR}/${SUMMARY_FILENAME} ${SUMMARY_DIR}/${SAMPLE_ID}.${SUMMARY_FILENAME}
done