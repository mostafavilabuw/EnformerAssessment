# Get ISM results for Enformer on a given gene set 
# For this analysis, get Enformer output for inserting each nucleotide at every position within 1000bp of the TSS

import time
import pandas as pd
import numpy as np
import sparse
import os
import tensorflow as tf
import tensorflow_hub as hub
import time as time
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

win = 1000 # window around TSS to get results for 

# enformer helper functs 
class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

def get_tss_ism(gene_id): 
    
    gene_res = np.zeros([win*2+1,4]) 
    
    model = Enformer(model_path)
    snp_info = pd.read_csv('/data/aspiro17/enformer_res/variant_info_100k/' + gene_id + '.csv',header=None,encoding='latin1')
    current_gene = np.load(save_path+'ref_seqs/'+gene_id+'.npy')

    starting_seq_len = np.shape(current_gene)[1]
    current_tss = int(gene_win_info['tss_hg38'][gene_id])
    attrib_start_pos = current_tss - int(input_len/2) + 1  

    # adjust sequence to be input to model
    current_gene = np.transpose(current_gene, (1,0)) # transpose to be seq_len x 4 
    current_gene = current_gene[:, [0,3,2,1]] # go from previous: seqP=='A',seqP=='T',seqP=='G',seqP=='C', to Enformer: 'ACGT'
    current_gene = np.reshape(current_gene, (1, starting_seq_len, 4)) # add a 1 dimen
    current_gene = np.pad(current_gene, pad_width=((0,0),(pad_before-start_index, pad_before-(starting_seq_len - end_index)), (0,0))) # pad seq 

    # from ref 
    ref_out = model.predict_on_batch(current_gene)
    ref_out = ref_out['human'][0]
    ref_val = np.sum(ref_out[center_bins,track_idx]) # 448 is center bin 

    np.save(save_path + 'attrib_res/TSS_win_ISM_res/' +gene_id  + '_ref_pred', ref_val)

    center_idx = (current_tss - attrib_start_pos) + pad_before # where in this array 
    window_start_pos = center_idx - win

    for i in range(win*2+1): 

        # fill in for ref 
        current_onehot = list(current_gene[:,window_start_pos+i,:][0]) 
        current_nuc_idx = nuc_onehots.index(current_onehot)
        gene_res[i,current_nuc_idx] = ref_val # this nuc represents no change from the ref 
        rel_nuc_idxs = list(range(4))
        rel_nuc_idxs.remove(current_nuc_idx) # don't have to re run for ref nuc

        for idx in rel_nuc_idxs: # the other 3 nucs
            inserted_seq = current_gene.copy()
            current_onehot = nuc_onehots[idx]
            inserted_seq[:,window_start_pos+i,:] = current_onehot # get current seq

            var_out = model.predict_on_batch(inserted_seq)
            var_out = var_out['human'][0]
            var_val = np.sum(var_out[center_bins,track_idx])  

            gene_res[i,idx] = var_val

        np.save(save_path + 'attrib_res/TSS_win_ISM_res/' + gene_id + '_' + str(win) + 'bp', gene_res)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('which_gpu', action="store", default='0')
    parser.add_argument('gene_file', action="store", default='test')

    args = parser.parse_args()
    
    os.environ["CUDA_VISIBLE_DEVICES"]=args.which_gpu
    genes = np.load(args.gene_file)
    
    for gene_id in genes: 
        print(gene_id)
        get_tss_ism(gene_id)


        
        
