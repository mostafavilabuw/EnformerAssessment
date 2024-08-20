This directory contains scripts used to process the genomic data used to make and analyze predictions for differential gene expression between 839 individuals using the Enformer model.

## Genotype_AMPAD_WGS_04_phasing_Eagle_script.sh
- create phased VCF files

## save_ref_seqs.py
- save onehot-encoded reference genomic sequences within a window (here, 100k) around the TSS from fasta files 
- fasta files can be downloaded from UCSC genomc browser, for ex, rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz ./

## save_snp_info.py
- from vcf, save variant information by gene for each SNP within a window (here, 100k) of the TSS

## extractSeq.py
- using reference genomic fasta files (downloaded as described in save_ref_seqs.py) and variant information by gene (saved as described in save_snp_info.py), create personalized one-hot encodings, saved by gene as sparse .npz files using the package sparse (https://anaconda.org/conda-forge/sparse)
- variant information by gene is accessed in snp = pd.read_csv(filepath+'variantNucleotide'+winSiz+'/'+data.iloc[j,0]+'.csv',header=None,encoding='latin1')
- when these one hot encodings are loaded using sparse (as in EnformerAssessment/enformer_analysis/basic_pred_gene_expr.py), they are of shape [8,200001,1161], where the first dimension of size 8 contains the 4 bases for paternal and the 4 bases for maternal sequences, the second dimension is the length, and the third 1161 individuals

