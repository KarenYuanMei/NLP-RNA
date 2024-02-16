#!/usr/bin/env Rscript

#Goal: to format Rsem output files into the correct format as input for EdgeR for differential expression analysis:

options(stringsAsFactors=F)

#set the working directory with the RSEM output files
#setwd("/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_052323_newseq_analysis/rsem_quant")

# Get the working directory from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  stop("No working directory specified", call. = FALSE)
}
working_dir <- args[1]
setwd(working_dir)

print (working_dir)
#find the files with the specified string pattern:
files = dir(pattern="*.RSEM_Quant.rsem.isoforms.results",recursive=T)

#read the rsem output file
inFile=read.delim(files[[1]],row.names=1)
genes = row.names(inFile)
counts=inFile$expected_count
tpm=inFile$TPM

for(i in 2:length(files)) {
  print(i)
  inFile=read.delim(files[[i]],row.names=1)
  inFile = inFile[match(genes,rownames(inFile)),]
  counts = cbind(counts, inFile$expected_count)
  tpm = cbind(tpm, inFile$TPM)
}


rownames(counts) = rownames(tpm) = substr(genes,1,18)
colnames(counts) = colnames(tpm) = gsub("/","_", gsub(".RSEM_Quant.isoforms.results","",files))

save(file="./RSEM_Quant.isoforms.counts.RData",counts)
save(file="./RSEM_Quant.isoforms.tpm.RData",tpm)