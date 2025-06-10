#!/bin/bash

if [[ -z "$1" ]]; then
    echo "Usage: sh run-jobs.sh [vdjpipe|igblast]"
    exit 1
fi

exclude="! -name '*.airr.json' ! -name '*test.json'"

if  [[ $1 == "vdjpipe" ]]; then
    pattern="job-vdjpipe*.json"
elif [[ $1 == "igblast" ]]; then
    pattern="job-igblast*.json"
else
    echo "Unknown mode: $1"
    echo "Valid modes are: vdjpipe, igblast"
    echo "Usage: sh run-jobs.sh [vdjpipe|igblast]"
    exit 1
fi

find . -path ./old-jobs -prune -o \( -type f -name "$pattern" $exclude \) -print | while read -r file; do
    echo "Submitting: $file"
    vdjserver-tools jobs submit "$file"
done