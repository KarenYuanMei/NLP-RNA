#!/bin/bash

#SBATCH --job-name=format_rsem_output_for_edgeR
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/analysis_code/
#SBATCH --output=%x.%j.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=5G
#SBATCH --time=2-04:00:00


#Goal: to format Rsem output files into the correct format as input for EdgeR for differential expression analysis:

#Rscript lithium_format_rsem_output.R
Rscript format_rsem_genes_output.R  

Rscript format_rsem_isoforms_output.R  

