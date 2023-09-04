This directory contains scripts used to process the genomic data used to make and analyze predictions for differential gene expression between 839 individuals using the Enformer model.


## save_ref_seqs.py
- save onehot-encoded reference genomic sequences within a window (here, 100k) around the TSS from fasta files 
- fasta files can be downloaded from UCSC genomc browser, for ex, rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz ./

## save_snp_info.py
- from vcf, save variant information by gene for each SNP within a window (here, 100k) of the TSS

## extractSeq.py
- using reference genomic fasta files (downloaded as described in save_ref_seqs.py) and variant information by gene (saved as described in save_snp_info.py), create personalized one-hot encodings, saved by gene as .npz files
