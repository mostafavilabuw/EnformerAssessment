# get ISM within 20bp of each driver SNP
# start with the most common variant inserted, then insert each possible variant
# same as from_main_var_drivers_ISM but with main variant inserted 

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

driver_win = 20
os.environ["CUDA_VISIBLE_DEVICES"]='0'

# enformer helper functs 
class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

def get_gene_drivers(gene_id, chrom): 
    
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
    
    # get all drivers pos
    driver_pos = driver_dict[gene_id]
    
    for snp in driver_pos: 

        snp = int(snp)
        driver_res = np.zeros([41,4]) # fill in the middle w ref
        ref = snp_info[snp_info[2] == snp][1].values[0] # get the single nuc

        # get most common alt 
        sel_row = snp_info[snp_info[2] == snp].iloc[:,3:] # 3 on 
        sel_row = sel_row.astype(str).values.flatten().tolist()
        sel_row = ' '.join(sel_row)
        counts = [sel_row.count('A'), sel_row.count('C'), sel_row.count('G'), sel_row.count('T')] 
        counts[nuc_order.index(ref)] = 0 # we want the most common allele that's not the ref 
        most_common_alt_idx = np.argmax(counts) # now take the max 
        most_common_alt_onehot = nuc_onehots[most_common_alt_idx]

        snp_center_idx = (snp - attrib_start_pos) + pad_before # where in this array 
        window_start_pos = snp_center_idx - driver_win
        
        current_gene_main_var_inserted = current_gene.copy()
        current_gene_main_var_inserted[:,snp_center_idx,:] = most_common_alt_onehot
        
        # from main var 
        main_var_out = model.predict_on_batch(current_gene_main_var_inserted)
        main_var_out = main_var_out['human'][0]
        main_var_val = np.sum(main_var_out[center_bins,track_idx]) bin 

        np.save(save_path + 'attrib_res/drivers_analysis/from_main_var/' +gene_id  + '_' + str(snp) + '_main_var_pred', main_var_val)
        
        for i in range(driver_win*2+1): 
            current_onehot = list(current_gene_main_var_inserted[:,window_start_pos+i,:][0])
            
            current_nuc_idx = nuc_onehots.index(current_onehot)
            driver_res[i,current_nuc_idx] = main_var_val # this nuc represents no change from the ref 
            rel_nuc_idxs = list(range(4))
            rel_nuc_idxs.remove(current_nuc_idx)
            
            for idx in rel_nuc_idxs: # the other 3
                inserted_seq = current_gene_main_var_inserted.copy()
                current_onehot = nuc_onehots[idx]                
                inserted_seq[:,window_start_pos+i,:] = current_onehot # get current seq
        
                var_out = model.predict_on_batch(inserted_seq)
                var_out = var_out['human'][0]
                var_val = np.sum(var_out[center_bins,track_idx]) # use 447,448,and 449 
                
                driver_res[i,idx] = var_val
        
        np.save(save_path + 'attrib_res/drivers_analysis/from_main_var/' +gene_id  + '_' + str(snp), driver_res)


if __name__ == '__main__':
    
    drivers = np.load('genes_to_run/drivers_to_run.npy')
    
    driver_dict = {}
    for item in drivers: 
        if item[0] not in driver_dict.keys(): 
            driver_dict[item[0]] = [item[1]]
        else:
            driver_dict[item[0]].append(item[1])
    
    for gene in driver_dict.keys():
        print(gene)
        chrom = all_gene_info[all_gene_info['gene_id'] == gene]['chr_hg38']
        get_gene_drivers(gene, chrom)