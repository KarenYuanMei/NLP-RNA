#Parse fasta file; analyze the sequences

#activate environment on the cluster: nlp-rna

import csv
import numpy as np
import pandas as pd

from collections import Counter

from statistics import mode, mean

from Bio import SeqIO

import pylab

df=pd.read_csv('~/Documents/nlp_mRNA/source_data/Ensembl_human_cdna/Ensembl_human_cdna.txt')
# print (df)
# print (df.head(50))


#annoying detail: the data source file has to be in the same folder as the analysis code:
# for seq_record in SeqIO.parse("Ensembl_human_cdna.txt", "fasta"):
# 	print(seq_record.id)
# 	print (seq_record.name)
# 	print(seq_record.seq)
# 	print(len(seq_record))

records = list(SeqIO.parse("Ensembl_human_cdna.txt", "fasta"))
print (len(records))

sizes = [len(rec) for rec in SeqIO.parse("Ensembl_human_cdna.txt", "fasta")]
#274081 sequences; 

print ('mode', mode(sizes))
print ('mean', mean(sizes))


data = Counter(sizes)
#print (data)

pylab.hist(sizes, bins=20)
pylab.title(
"%i human cdna sequences\nLengths %i to %i" % (len(sizes), min(sizes), max(sizes))
)
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.yscale("log")  


image_format = 'svg' # e.g .png, .svg, etc.
image_name = 'myimage.svg'
pylab.savefig(image_name, format=image_format)