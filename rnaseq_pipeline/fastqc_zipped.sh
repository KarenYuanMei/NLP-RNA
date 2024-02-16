#!/bin/bash

#SBATCH --job-name=fastqc_array
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --output=%x.%j.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=2-04:00:00
#SBATCH --array=1-${NUM_FILES}

# Accept FASTQ_DIR and NUM_FILES from arguments
FASTQ_DIR=$1
NUM_FILES=$2
OUTPUT_DIR=$3


# Get all fastq.gz files
files=("$FASTQ_DIR"/*.fastq.gz)

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"


# Use the array task ID to get the file to process
FILE="${files[$SLURM_ARRAY_TASK_ID-1]}"


# Define the output file path
output_file="${OUTPUT_DIR}/$(basename "$FILE" .fastq.gz)_fastqc.zip"

# Run FastQC if the output file doesn't exist
if [[ ! -f $output_file ]]; then
    fastqc "$FILE" -o "$OUTPUT_DIR"
fi

