# Test affect of data augmentation (shifts and reverse complement) on Enformer output

import time
import pandas as pd
import numpy as np
import sparse
import os

import tensorflow as tf
import tensorflow_hub as hub

import time as time

import kipoiseq
from kipoiseq import Interval

save_path = '/data/aspiro17/enformer_res/'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
data_path = '/data/mostafavilab/bng/rosmapAD/projects/insilicoMutagenesis/extractSequence/results/sequence100K/'

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

shifts = [-3, -2, -1, 0, 1, 2, 3] # how to shift sequence

# enformer helper functs 
class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}


# other helper functs

def pad_genes(pat_current_gene, mat_current_gene,shift=0):
    # pad genes to create a shift 

    padded_pat_current_gene = np.pad(pat_current_gene, pad_width=((0,0),(pad_before-start_index+shift, pad_before-(starting_seq_len - end_index)-shift), (0,0)))
    padded_mat_current_gene = np.pad(mat_current_gene, pad_width=((0,0),(pad_before-start_index+shift, pad_before-(starting_seq_len - end_index)-shift), (0,0)))
    
    return padded_pat_current_gene, padded_mat_current_gene


def rev_comp(pat_current_gene, mat_current_gene):
    
    pat_comp = pat_current_gene[:,:,[3,2,1,0]]
    mat_comp = mat_current_gene[:,:,[3,2,1,0]]
    
    rc_pat = np.flip(pat_comp, axis=1)
    rc_mat = np.flip(mat_comp, axis=1)
    
    return rc_pat, rc_mat


def get_aug_pred(gene_id, curr_chr): 
    
    # data for specific subjects 
    curr_gene_path = data_path + 'chr' + str(curr_chr) + '/' + gene_id + '.npz'
    current_sparse = sparse.load_npz(curr_gene_path)
    current_gene = current_sparse.todense()
    
    # if not running for all subj
    current_gene = np.take(current_gene, final_idx, axis=2)
    num_subj = current_gene.shape[2]
    
    # init results 
    gene_res = np.zeros([num_subj,2,8],dtype=float)

    pat_current_gene = current_gene[:4,:]
    mat_current_gene = current_gene[4:,:]
    
    # transpose to be num_sub x seq_len x 4 
    pat_current_gene = np.transpose(pat_current_gene, (2,1,0))
    mat_current_gene = np.transpose(mat_current_gene, (2,1,0))

    # go from previous: seqP=='A',seqP=='T',seqP=='G',seqP=='C', to Enformer: 'ACGT'    pat_current_gene = pat_current_gene[:, :, [0,3,2,1]]
    mat_current_gene = mat_current_gene[:, :, [0,3,2,1]]
    
    # get output for shifted sequence 
    for shift in range(len(shifts)): 
        padded_pat_current_gene, padded_mat_current_gene = pad_genes(pat_current_gene, mat_current_gene,shift=shift)
        
        for sub in range(num_subj):
            
            pat_single_sub = np.reshape(padded_pat_current_gene[sub,:,:], (1, padded_input_len, 4))
            mat_single_sub = np.reshape(padded_mat_current_gene[sub,:,:], (1, padded_input_len, 4))

            pat_out = model.predict_on_batch(pat_single_sub)
            mat_out = model.predict_on_batch(mat_single_sub)
                
            pat_out = pat_out['human'][0]
            mat_out = mat_out['human'][0]
            
            pat_bins_sum = pat_out[center_bins[0]][track_idx] + pat_out[center_bins[1]][track_idx] + pat_out[center_bins[2]][track_idx]
            mat_bins_sum = mat_out[center_bins[0]][track_idx] + mat_out[center_bins[1]][track_idx] + mat_out[center_bins[2]][track_idx]
            
            gene_res[sub,0,shift] = pat_bins_sum
            gene_res[sub,1,shift] = mat_bins_sum
    
    # get output for rc sequence 
    rc_pat, rc_mat = rev_comp(pat_current_gene, mat_current_gene)
    padded_pat_current_gene, padded_mat_current_gene = pad_genes(rc_pat, rc_mat,shift=0)

    for sub in range(num_subj):
        pat_single_sub = np.reshape(padded_pat_current_gene[sub,:,:], (1, 393216, 4))
        mat_single_sub = np.reshape(padded_mat_current_gene[sub,:,:], (1, 393216, 4))

        pat_out = model.predict_on_batch(pat_single_sub)
        mat_out = model.predict_on_batch(mat_single_sub)

        pat_out = pat_out['human'][0]
        mat_out = mat_out['human'][0]

        pat_bins_sum = pat_out[center_bins[0]][track_idx] + pat_out[center_bins[1]][track_idx] + pat_out[center_bins[2]][track_idx]
        mat_bins_sum = mat_out[center_bins[0]][track_idx] + mat_out[center_bins[1]][track_idx] + mat_out[center_bins[2]][track_idx]

        gene_res[sub,0,7] = pat_bins_sum
        gene_res[sub,1,7] = mat_bins_sum
            
    np.save(save_path + 'data_aug_pred/' + gene_id, gene_res)

            
if __name__ == '__main__':
    
    # deal with subject indexing 
    
    # all subj 
    chosen_subs = np.load('/data/aspiro17/enformer_res/sub_lists/all_exp_and_seq_subs.npy')
    
    # select a random 100 subjs
    #chosen_subs = np.load(save_path+ 'sub_lists/second_100_subs.npy')
    
    vcf_sub_list = np.load(save_path + 'sub_lists/vcf_subs.npy')

    # get indices for randomly chosen subs within the 1161 subs we have vcf data for 
    idx_in_vcf = []
    for item in chosen_subs: 
        idx_in_vcf.append(np.where(vcf_sub_list == item))

    final_idx = []
    for item in idx_in_vcf: 
        final_idx.append(item[0][0])
        
    # load model
    model = Enformer(model_path)

    # load TSS info
    gene_win_info = pd.read_csv('/data/mostafavilab/bng/rosmapAD/data/expressionData/DLPFC/20220207-bulk-RNAseq/gene-ids-and-positions.tsv', sep = '\t', index_col = 1)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('which_gpu', action="store", default='0')
    parser.add_argument('gene_file', action="store", default='test')
    args = parser.parse_args()
    
    os.environ["CUDA_VISIBLE_DEVICES"]=args.which_gpu
    genes = np.load(args.gene_file)
    
    for gene in genes:
        print(gene)
        chrom = gene_win_info[gene_win_info['gene_id'] == gene]['chr_hg38']
        get_aug_pred(gene, chrom)
    
    
    