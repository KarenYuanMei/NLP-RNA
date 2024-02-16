#!/bin/bash

#SBATCH --job-name=rsem_quant
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --output=%x.%A.%a.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=4-04:00:00
#SBATCH --array=1-${PAIRED_FILES}  

# Accept TRIMMED_DIR and NUM_FILES from arguments
ALIGNED_DIR=$1
PAIRED_FILES=$2
RSEM_OUTPUT_DIR=$3


#for a batch of files:

#paths to required files:
#path to trimmed fastqs:


# Create the output directory if it doesn't exist
mkdir -p "$RSEM_OUTPUT_DIR"

#path to RSEM index:
rsem_index=/cellar/users/yumei/Documents/mouse_genome_assembly/gencode_genomes/rsem_ref/mouse_gencode

#set the path to RSEM

export PATH=/cellar/users/yumei/Documents/RSEM-1.3.3/bin:$PATH



# Get a list of files
files=($ALIGNED_DIR/*out_sorted.bam*)

input_file="${files[$SLURM_ARRAY_TASK_ID-1]}"

echo $input_file

base_file=$(basename $input_file "_Aligned.out_sorted.bam")
new_base=${base_file#"star_aligned_files/"}

#prefix="${filename%%_S*}"
#output_prefix=${new_base:0:3}
output_prefix="${new_base%%_Aligned*}"


echo $output_prefix

rsem-calculate-expression  -p 12 --paired-end --bam --estimate-rspd --calc-ci --no-bam-output --ci-memory 1000 --seed 12345 --forward-prob 0 $input_file  $rsem_index  $RSEM_OUTPUT_DIR/${output_prefix}_rsem.RSEM_Quant.rsem >& $RSEM_OUTPUT_DIR/Log_out.rsem

# done



