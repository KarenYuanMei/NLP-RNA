#!/bin/bash

#SBATCH --job-name=mmseqs2
#SBATCH --chdir=/cellar/users/yumei/Documents/nlp_mRNA/P2_RNAseq_finetuning_by_gene/
#SBATCH --output=%x.%j.out
#SBATCH --partition=nrnb-compute
#SBATCH --account=nrnb-compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=10-09:00:00 



# Convert both the training and validation FASTA files into MMseqs2 databases:
mmseqs createdb P2_cortex_finetuning_by_gene_nontest.fasta nontestDB
mmseqs createdb P2_cortex_finetuning_by_gene_test.fasta testDB

tmpDir="/cellar/users/yumei/Documents/nlp_mRNA/P2_RNAseq_finetuning_by_gene/tmpDir"

# Use MMseqs2 to search the sequences from the training database against the validation database:
mmseqs search nontestDB testDB resultDB $tmpDir --min-seq-id 0.5 --alignment-mode 3 --max-seqs 300 -s 7 -c 0.8 --cov-mode 0 --search-type 3

# # Convert the MMseqs2 results into a human-readable format:
mmseqs convertalis nontestDB testDB resultDB result.tsv

# # Clean up
rm -rf $tmpDir

awk '$3 >= 0.5 {print $1}' result.tsv > sequences_to_remove.txt

sort sequences_to_remove.txt | uniq > unique_sequences_to_remove.txt


seqkit grep -v -f unique_sequences_to_remove.txt P2_cortex_finetuning_by_gene_nontest.fasta > P2_cortex_by_gene_filtered_nontest.fasta
