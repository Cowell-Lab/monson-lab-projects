#!/bin/bash

find . -type f -name "*.json" ! -name "*.airr.json" ! -name "*test.json" | while read -r file; do
    echo "Submitting: $file"
    vdjserver-tools jobs submit $file
done
