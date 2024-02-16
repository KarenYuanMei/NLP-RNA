#!/bin/bash

#SBATCH --job-name=cutadapt_trim
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --output=%x.%j.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=2-04:00:00
#SBATCH --array=1-${PAIRED_FILES}

# Accept TRIMMED_DIR and NUM_FILES from arguments
FASTQ_DIR=$1
PAIRED_FILES=$2
TRIMMED_DIR=$3


#Goal: to use Cutadapt 3.5 to trim away any adaptor sequences and low quality reads

#activate system
#module load slurm

# Create the output directory if it doesn't exist
mkdir -p "$TRIMMED_DIR"


# Get all prefixes
prefixes=()
for file in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
  if [ -f "$file" ]; then
    filename=$(basename "$file")
    prefix=$(echo "$filename" | sed 's/_merged_R1_001.fastq.gz//')
    if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]; then
        prefixes+=("$prefix")
    fi
  fi
done

# Get the prefix for this task
prefix="${prefixes[$SLURM_ARRAY_TASK_ID-1]}"

# Find R1 and R2 files for the current prefix
R1_files=$(find "$FASTQ_DIR" -type f -name "${prefix}_merged_R1_001.fastq.gz")
R2_files=$(find "$FASTQ_DIR" -type f -name "${prefix}_merged_R2_001.fastq.gz")

echo $prefix
echo "For prefix $prefix, found R1 files: $R1_files"
echo "For prefix $prefix, found R2 files: $R2_files"


for R1 in $R1_files; do
  R2="${R1/R1/R2}"
  if [[ -e $R2 ]]; then
    echo "Processing $R1 and $R2 with cutadapt..."
    cutadapt --minimum-length=5 -q 30 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $TRIMMED_DIR/$(basename $R1 .fastq.gz).trim.fastq.gz -p $TRIMMED_DIR/$(basename $R2 .fastq.gz).trim.fastq.gz $R1 $R2 
  else
    echo "Matching file for $R1 not found"
  fi
done


