import numpy as np
import sys, os

gene = sys.argv[1]
snpfile = np.load(gene+'snp_info.npz', allow_pickle = True)
snps = snpfile['snps']
locs = snpfile['rows']

freq = np.around(np.mean(snps>0, axis = 1),3)
if '--centered' in sys.argv:
    freq[freq>0.5] = np.around(1.-freq[freq>0.5],3)
    np.savetxt(gene+'_frequencycentered.txt', np.array([locs, freq]).T, fmt = '%s')
else:
    np.savetxt(gene+'_frequency.txt', np.array([locs, freq]).T, fmt = '%s')


