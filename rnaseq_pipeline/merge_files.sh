#!/bin/bash

# Accept TRIMMED_DIR and NUM_FILES from arguments
FASTQ_singles=$1
FASTQ_DIR=$2

# Get all prefixes
prefixes=()
for file in "$FASTQ_singles"/*; do
  if [ -f "$file" ]; then
    filename=$(basename "$file")
    prefix=$(echo "$filename" | sed 's/_S[0-9]\+.*//')
    # Append prefix to array if it does not exist
    if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]; then
        prefixes+=("$prefix")
    fi
  fi
done

for prefix in "${prefixes[@]}"; do
    echo "Processing prefix: $prefix"
    
    # Find and merge R1 files
    cat "$FASTQ_singles"/${prefix}*_R1_001.fastq.gz > "$FASTQ_DIR"/${prefix}_merged_R1_001.fastq.gz
    
    # Find and merge R2 files
    cat "$FASTQ_singles"/${prefix}*_R2_001.fastq.gz > "$FASTQ_DIR"/${prefix}_merged_R2_001.fastq.gz
done

