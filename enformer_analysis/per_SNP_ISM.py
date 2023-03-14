# Get ISM results for Enformer on a given gene set 
# For this analysis, get Enformer output for reference and Enformer output for inserting the most common variant for each SNP position in a gene 

import time
import pandas as pd
import numpy as np
import sparse
import os
import tensorflow as tf
import tensorflow_hub as hub
import argparse

track_idx = 4980 # CAGE:brain, adult,
center_bins = [447,448,449] # three center bins of enformer output

save_path = '/data/aspiro17/enformer_res/'
model_path = 'https://tfhub.dev/deepmind/enformer/1'

padded_input_len = 393216 # full input to enformer 
input_len = 196608
starting_seq_len = 200001 # starting length of reference seq 
pad_before = int((padded_input_len - input_len)/2) # "pad" sequence to padded_input_len

mid_index = int((starting_seq_len-1)/2)
start_index = int(mid_index - (input_len-1)/2)
end_index = int(mid_index + (input_len-1)/2) + 1 

nuc_order =['A', 'C', 'G', 'T']
nuc_onehots = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

gene_win_info = pd.read_csv('/data/mostafavilab/bng/rosmapAD/data/expressionData/DLPFC/20220207-bulk-RNAseq/gene-ids-and-positions.tsv', sep = '\t', index_col = 1) # contains TSS position info 


# enformer helper functs 
class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

def run_ism(gene_id, chrom): 
    # for a current gene_id (ENSG) and chromosome, save ISM results 
    model = Enformer(model_path)
    
    snp_info = pd.read_csv('/data/aspiro17/enformer_res/variant_info_100k/' + gene_id + '.csv',header=None,encoding='latin1') # variant info for this gene 
    current_gene = np.load(save_path+'ref_seqs/'+gene_id+'.npy') # reference sequence for this gene 
    
    most_common_alts = [] # save most common alternate alleles for future analysis 
    var_vals = [] # save output of enformer for inserting a given var (for each var position) 

    starting_seq_len = np.shape(current_gene)[1]
    current_tss = int(gene_win_info['tss_hg38'][gene_id])
    attrib_start_pos = current_tss - int(input_len/2) + 1  

    # get snp positions for this gene (within the 100k window) 
    all_snp_pos = list(snp_info[2])
    snp_pos = []
    for snp in all_snp_pos: 
        adjusted_pos = snp - attrib_start_pos # make the attrib start pos like 0 
        if adjusted_pos >= 0 and adjusted_pos <= input_len: # in range
            snp_pos.append(snp)
    np.save(save_path + 'snp_positions/' + gene_id + '_snp_positions', snp_pos)

    # adjust sequence to be input to model
    current_gene = np.transpose(current_gene, (1,0)) # transpose to be seq_len x 4 
    current_gene = current_gene[:, [0,3,2,1]] # go from previous: seqP=='A',seqP=='T',seqP=='G',seqP=='C', to Enformer: 'ACGT'
    current_gene = np.reshape(current_gene, (1, starting_seq_len, 4)) # add a 1 dimen
    current_gene = np.pad(current_gene, pad_width=((0,0),(pad_before-start_index, pad_before-(starting_seq_len - end_index)), (0,0))) # pad seq 
    
    ref_out = model.predict_on_batch(current_gene)
    ref_out = ref_out['human'][0]
    ref_val = np.sum(ref_out[center_bins,track_idx]) # sum over center 3 bins 
    
    np.save(save_path + 'attrib_res/ism_res/' +gene_id  + '_ref_pred', ref_val)
    
    for snp in snp_pos: 
        print(snp)
        
        inserted_seq = current_gene.copy()
        ref = snp_info[snp_info[2] == snp][1].values[0] # single ref nuc 
        
        # get most common alt 
        sel_row = snp_info[snp_info[2] == snp].iloc[:,3:] # 3 on 
        sel_row = sel_row.astype(str).values.flatten().tolist()
        sel_row = ' '.join(sel_row)
        counts = [sel_row.count('A'), sel_row.count('C'), sel_row.count('G'), sel_row.count('T')] 
        counts[nuc_order.index(ref)] = 0 # we want the most common allele that's not the ref 
        alt_idx = np.argmax(counts) # now take the max 
        allele_to_insert = nuc_onehots[alt_idx]
        
        most_common_alt = nuc_order[alt_idx]
        most_common_alts.append(most_common_alt)

        idx_to_insert = (snp - attrib_start_pos) + pad_before # where to insert variant 
        
        print(ref, current_gene[0,idx_to_insert,:]) # make sure the reference value is the same as what we find at this idx in the one-hot
        
        inserted_seq[:,idx_to_insert,:] = allele_to_insert # insert 
        
        var_out = model.predict_on_batch(inserted_seq)
        var_out = var_out['human'][0]
        var_val = np.sum(var_out[center_bins,track_idx]) 

        var_vals.append(var_val)                     
            
    np.save(save_path + 'attrib_res/ism_res/' +gene_id  + '_var_preds', var_vals)
    np.save(save_path + 'attrib_res/ism_res/' +gene_id  + '_alt_alleles', most_common_alts)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('which_gpu', action="store", default='0')
    parser.add_argument('gene_file', action="store", default='test')
    args = parser.parse_args()
    
    os.environ["CUDA_VISIBLE_DEVICES"]=args.which_gpu
    genes = np.load(args.gene_file)
    
    for gene in genes:
        print(gene)
        chrom = gene_win_info[gene_win_info['gene_id'] == gene]['chr_hg38']
        run_ism(gene, chrom)