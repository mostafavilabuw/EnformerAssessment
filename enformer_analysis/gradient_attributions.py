# Save Enformer gradients 
# Save (gradient at reference sequence x reference) and (gradient at reference sequence x main variant) for each main variant (most common alternate allele) at each SNP position for a given gene

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
# enformer helper functs 
class Enformer:

  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass
    grad = tape.gradient(prediction, input_sequence)
    return grad


def get_ref_attribs(gene_id, chrom): 
    
    model = Enformer(model_path)
    
    snp_info = pd.read_csv('/data/aspiro17/enformer_res/variant_info_100k/' + gene_id + '.csv',header=None,encoding='latin1')

    most_common_alts = []
    current_gene = np.load(save_path+'ref_seqs/'+gene_id+'.npy')    
    
    starting_seq_len = np.shape(current_gene)[1]
    current_tss = int(gene_win_info['tss_hg38'][gene_id])
    attrib_start_pos = current_tss - int(input_len/2) + 1  
    
    snp_pos = np.load(save_path + 'snp_positions/' + gene_id + '_snp_positions.npy')

    # adjust sequence to be input to model
    current_gene = np.transpose(current_gene, (1,0)) # transpose to be seq_len x 4 
    current_gene = current_gene[:, [0,3,2,1]] # go from previous: seqP=='A',seqP=='T',seqP=='G',seqP=='C', to Enformer: 'ACGT'
    current_gene = np.reshape(current_gene, (1, starting_seq_len, 4)) # add a 1 dimen
    current_gene = np.pad(current_gene, pad_width=((0,0),(pad_before-start_index, pad_before-(starting_seq_len - end_index)), (0,0))) # pad seq 

    inserted_seq = current_gene.copy()
    
    for snp in snp_pos: 
        
        ref = snp_info[snp_info[2] == snp][1].values[0] 

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
        
        idx_to_insert = (snp - attrib_start_pos) + pad_before
        inserted_seq[:,idx_to_insert,:] = allele_to_insert 
        
    np.save(save_path + 'most_common_alt_alleles/' +gene_id  + '_alt_alleles', most_common_alts)  
        
    # target mask to get relevant gradient 
    target_mask = np.zeros([896,5313], dtype='float32')
    for idx in [447, 448, 449]: # central 3 bins 
        target_mask[idx, track_idx] = 1

    # get ref grad 
    ref_grad = model.contribution_input_grad(current_gene[0,:,:].astype(np.float32), target_mask).numpy() 
    np.save(save_path + 'attrib_res/ref_attribs/' + gene_id + '_complete_grad_at_ref', ref_grad)

    # get (grad at ref) x (ref) 
    grad_at_ref_times_ref = ref_grad * current_gene[0,:,:]
    grad_at_ref_times_ref = tf.squeeze(grad_at_ref_times_ref, axis=0)
    grad_at_ref_times_ref = tf.reduce_sum(grad_at_ref_times_ref, axis=-1).numpy()
    
    # get (grad at ref) x (var) 
    grad_at_ref_times_var = ref_grad * inserted_seq
    grad_at_ref_times_var = tf.squeeze(grad_at_ref_times_var, axis=0)
    grad_at_ref_times_var = tf.reduce_sum(grad_at_ref_times_var, axis=-1).numpy()

    grad_times_ref_vals = []
    grad_times_var_vals = []
    for snp in snp_pos: 
        adjusted_pos = (snp - attrib_start_pos) + pad_before # where in the array 
        grad_times_ref_vals.append(grad_at_ref_times_ref[adjusted_pos]) 
        grad_times_var_vals.append(grad_at_ref_times_var[adjusted_pos])

    np.save(save_path + 'attrib_res/ref_attribs/' + gene_id + '_grad_at_ref_times_ref', grad_times_ref_vals)
    np.save(save_path + 'attrib_res/ref_attribs/' + gene_id + '_grad_at_ref_times_var', grad_times_var_vals)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('which_gpu', action="store", default='0')
    parser.add_argument('gene_file', action="store", default='test')

    args = parser.parse_args()
    
    os.environ["CUDA_VISIBLE_DEVICES"]=args.which_gpu
    genes = np.load(args.gene_file)
    
    for gene in genes:
        print(gene)
        chrom = all_gene_info[all_gene_info['gene_id'] == gene]['chr_hg38']
        get_ref_attribs(gene, chrom)
    


