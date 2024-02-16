#!/bin/bash

#SBATCH --job-name=rna-seq_master
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --time=4-04:00:00


#check the format of the files: Does each R1 file have multiple files/lanes? Does each R2 file have multuiple files/lanes? If so, then should merge to form
#one R1 file and one R2 file for each sample


# Set the directory containing the .fastq.gz files
#FASTQ_singles='/cellar/users/yumei/Documents/nlp_mRNA/experimental_raw_data/mouse_lji_110623/mouse_lji_110623_singles'


FASTQ_DIR='/cellar/users/yumei/Documents/nlp_mRNA/experimental_raw_data/mouse_lji_110623/mouse_lji_110623_merged'


#sbatch merge_files.sh "${FASTQ_singles}" "${FASTQ_DIR}" 


#Step 1: Fastqc===================================================================================================

source /cellar/users/yumei/miniconda3/etc/profile.d/conda.sh

#environment: fastqc_env
conda activate fastqc_env

NUM_FILES=$(ls -1 "${FASTQ_DIR}"/*.fastq.gz | wc -l)

# OUTPUT_DIR="/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/fastqc_zipped_files"


# # Submit FastQC jobs and capture the job ID

# FASTQC_JOB_ID=$(sbatch --parsable --array=1-${NUM_FILES} fastqc_zipped.sh "${FASTQ_DIR}" "${NUM_FILES}" "${OUTPUT_DIR}")

# # Wait for the FastQC jobs to complete
# echo "Waiting for FastQC jobs to complete..."
# while squeue -j $FASTQC_JOB_ID | grep -q "$FASTQC_JOB_ID"; do 
#     sleep 60  # Wait for 60 seconds before checking again
# done

# echo "FastQC jobs completed."

# multiqc $OUTPUT_DIR/*


# conda deactivate


# #Step 2: trim the adapters using Cutadapt ===================================================================
# #environment: cutadaptenv
# conda activate cutadaptenv

PAIRED_FILES=$((NUM_FILES / 2))


# #path to output folder: trimmed fastqs:
TRIMMED_DIR='/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/trimmed_fastqs/'


# # Call the fastqc_zipped.sh script with arguments and set the SLURM array size
# #sbatch --array=1-${PAIRED_FILES} cutadapt.sh "${FASTQ_DIR}" "${PAIRED_FILES}" "${TRIMMED_DIR}"

# TRIM_JOB_ID=$(sbatch --parsable --array=1-${PAIRED_FILES} cutadapt.sh "${FASTQ_DIR}" "${PAIRED_FILES}" "${TRIMMED_DIR}")

# # Wait for the FastQC jobs to complete
# echo "Waiting for TRIM jobs to complete..."
# while squeue -j $TRIM_JOB_ID | grep -q "$TRIM_JOB_ID"; do 
#     sleep 600  # Wait for 60 seconds before checking again
# done

# echo "TRIM jobs completed."

# conda deactivate


#Step 3: fastqc the trimmed fastqs============================================================================

#conda activate fastqc_env

# Submit FastQC jobs and capture the job ID
#FASTQC_TRIMMED_DIR="/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/fastqc_trimmed_files"


# FASTQC_JOB_ID=$(sbatch --parsable --array=1-${NUM_FILES} fastqc_zipped.sh "${TRIMMED_DIR}" "${NUM_FILES}" "${FASTQC_TRIMMED_DIR}")

# # Wait for the FastQC jobs to complete
# echo "Waiting for FastQC jobs to complete..
# while squeue -j $FASTQC_JOB_ID | grep -q "$FASTQC_JOB_ID"; do 
#     sleep 60  # Wait for 60 seconds before checking again
# done

# echo "FastQC jobs completed."

# multiqc $FASTQC_TRIMMED_DIR/*

# conda deactivate

#Step 4: Star align the trimmed fastqs=====================================================================

conda activate cutadaptenv

ALIGNED_DIR="/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/star_aligned_files"


# #ALIGN_JOB_ID=$(sbatch --parsable --array=1-${PAIRED_FILES} star_align.sh "${TRIMMED_DIR}" "${PAIRED_FILES}" "${ALIGNED_DIR}")
ALIGN_JOB_ID=$(sbatch --parsable --array=1-${PAIRED_FILES} star_align.sh "${TRIMMED_DIR}" "${PAIRED_FILES}" "${ALIGNED_DIR}")

# Wait for the FastQC jobs to complete
echo "Waiting for ALIGN jobs to complete..."
while squeue -j $ALIGN_JOB_ID | grep -q "$ALIGN_JOB_ID"; do 
    sleep 600  # Wait for 60 seconds before checking again
done

echo "ALIGN jobs completed."


#Step 5: RSEM quantify read count ==========================================================================

RSEM_OUTPUT_DIR="/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/rsem_output"

RSEM_JOB_ID=$(sbatch --parsable --array=1-${PAIRED_FILES} rsem_quant.sh "${ALIGNED_DIR}" "${PAIRED_FILES}" "${RSEM_OUTPUT_DIR}")

# Wait for the FastQC jobs to complete
echo "Waiting for RSEM jobs to complete..."
while squeue -j $RSEM_JOB_ID | grep -q "$RSEM_JOB_ID"; do 
    sleep 600  # Wait for 60 seconds before checking again
done

echo "RSEM jobs completed."

conda deactivate

#Step 6: format rsem_output ===============================================================================


Rscript format_rsem_genes_output.R $RSEM_OUTPUT_DIR

Rscript format_rsem_isoforms_output.R $RSEM_OUTPUT_DIR


