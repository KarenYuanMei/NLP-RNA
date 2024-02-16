#!/bin/bash

# Accept TRIMMED_DIR from arguments
TRIMMED_DIR=$1
ALIGNED_DIR=$2

#SBATCH --job-name=star_trimmed
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --output=%x.%j.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=400G
#SBATCH --time=2-04:00:00

# Create the output directory if it doesn't exist
mkdir -p "$ALIGNED_DIR"

# Set samtools path
export PATH=/cellar/users/yumei/Documents/samtools-1.14/bin:$PATH


# # Get all prefixes
# prefixes=()
# for file in "$TRIMMED_DIR"/*_R1_001.trim.fastq.gz; do
#   if [ -f "$file" ]; then
#     filename=$(basename "$file")
#     prefix=$(echo "$filename" | sed 's/_merged_R1_001.trim.fastq.gz//')
#     if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]; then
#         prefixes+=("$prefix")
#     fi
#   fi
# done

# for prefix in "${prefixes[@]}"; do
# # Find R1 and R2 files for the current prefix
# R1_file=$(find "$TRIMMED_DIR" -type f -name "${prefix}_merged_R1_001.trim.fastq.gz")
# R2_file=$(find "$TRIMMED_DIR" -type f -name "${prefix}_merged_R2_001.trim.fastq.gz")
  
# echo $R1_file
# echo $R2_file
# echo $prefix

# ~/Documents/STAR/source/STAR --runThreadN 20 --genomeDir ~/Documents/mouse_genome_assembly/genome_index/ --readFilesIn $R1_file $R2_file  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within --outSAMattributes NH HI AS NM MD --outFilterMismatchNoverReadLmax 0.04 --sjdbScore 1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM  --outFileNamePrefix $ALIGNED_DIR/${prefix}_
# cat <( samtools view -H $ALIGNED_DIR/${prefix}_Aligned.toTranscriptome.out.bam ) <( samtools view -@ 1 $ALIGNED_DIR/${prefix}_Aligned.toTranscriptome.out.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S 10000000 -T ./ | tr ' ' '\n' ) | samtools view -@ 1 -bS - > $ALIGNED_DIR/${prefix}_Aligned.toTranscriptome.out_sorted.bam

# done


R1_file=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/trimmed_fastqs/Cytsol59A_merged_R1_001.trim.fastq.gz

R2_file=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/trimmed_fastqs/Cytsol59A_merged_R2_001.trim.fastq.gz

~/Documents/STAR/source/STAR --runThreadN 5 --genomeDir ~/Documents/mouse_genome_assembly/genome_index/ --readFilesIn $R1_file $R2_file  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within --outSAMattributes NH HI AS NM MD --outFilterMismatchNoverReadLmax 0.04 --sjdbScore 1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM  --outFileNamePrefix $ALIGNED_DIR/${prefix}_
