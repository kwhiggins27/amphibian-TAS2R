#!/bin/bash

int=$1
file=$2

# Read the file and perform sbatch command for the first $int items
count=0
while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ $count -ge $int ]]; then
        break
    fi

    sbatch "./pipeline_v7.sh" "$line"

    # sleep 60
    ((count++))

    sleep 1
done < "$file"
