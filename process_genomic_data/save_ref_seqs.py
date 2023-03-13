# For relevant genes, save reference sequence data from hg38 within 100k of of the TSS into one-hot encoded form 

import pandas as pd
import numpy as np
import gzip
import sparse
import os
import argparse

WIN = 1e5
NUM_SYM = 4

win_path='/data/mostafavilab/bng/rosmapAD/projects/insilicoMutagenesis/extractSequence/results/geneWin100K.txt'
save_path = '/data/aspiro17/enformer_res/ref_seqs/'
ref_path = '/data/aspiro17/seqtoexp/hg38/' # using data downloaded from UCSC genome browser 

gene_win_info = pd.read_csv('/data/mostafavilab/bng/rosmapAD/projects/insilicoMutagenesis/extractSequence/results/geneWin100K.txt', sep='\t', header=None)
all_gene_ids = np.loadtxt('/data/aspiro17/enformer_res/gene_lists/genes_to_run.txt',dtype=str) # genes to process 

for CHR in range(1,23):
    print(CHR)
    
    # load ref seq 
    with gzip.open(ref_path+'chr'+str(CHR)+'.fa.gz','rt') as f: 
        f.readline() # Remove non-sequence line
        ref = f.read().replace('\n','')
    seqR = np.array(list(ref.upper()))
    
    print('chrom loaded')
    
    chrom_len = len(seqR)
    
    # for relevant genes on this chrom 
    rel_genes = gene_win_info[gene_win_info[1] == CHR]
    num_genes = rel_genes.shape[0]
        
    for i in range(num_genes):
        print(i)
        current_gene = seqR[int(rel_genes.iloc[i][2]):int(rel_genes.iloc[i][3])+1] # take the relevant window 
        current_onehot = np.zeros((NUM_SYM,int(2*WIN+1)),dtype='i8')  
        current_onehot[:,:] = [current_gene=='A',current_gene=='T',current_gene=='G',current_gene=='C'] # note: this order is different from what is used in Enformer, we adjust the onehot later to match that order 

        np.save(save_path+rel_genes.iloc[i,0]+'.npy',current_onehot)
        
    
