# Uses mean and std of absolute correlation from Null and compute p-value and corrected p-value for each gene
# Usage
# python3 compute_tstatpvalue.py ../Prediction_correlationsCageAdultBrain.txt GeneSpecific_CorrelationtoObsRandomNull.txt

import numpy as np
import scipy.stats as stats
import sys, os
from statsmodels.stats.multitest import multipletests
# Read file with correlation between observed and predicted expression for all individuals
corr = np.genfromtxt(sys.argv[1], dtype = str)
# Reat file with mean and std of random null
null = np.genfromtxt(sys.argv[2], dtype = str)

# sort corr file so that genes match null file
cnames, csort = np.unique(corr[:,0], return_index = True)
csort = csort[np.isin(cnames, null[:,0])]
corr = corr[csort]
# sort null file so that genes match corr file
nsort = np.argsort(null[:,0])[np.isin(np.sort(null[:,0]), corr[:,0])]
null =null[nsort]

# check sorting process
print(np.array_equal(corr[:,0], null[:,0]))

# select columns that contain the data: 
# sm: predicted correlation
# m: mean of Null
# sv: std of Null
# n: number of random sets to compute mean and std
sm, m, sv, n = np.absolute(corr[:,1].astype(float)), null[:,1].astype(float), null[:,2].astype(float), null[:,3].astype(float)
sm, m, sv = np.nan_to_num(sm), np.nan_to_num(m), np.nan_to_num(sv)

tt = (sm-m)/np.sqrt(sv/n)  # t-statistic for mean
pval = stats.t.sf(np.abs(tt), n-1)
if  '--bothsided' in sys.argv:
    pval = pval *2
else:
    pval[tt < 0] = 1

# Correct for multiple testing with benjamini hochberg
issig, corr_pvals, calphasid, calphabonf = multipletests(np.nan_to_num(pval,nan=1), alpha=0.05, method='fdr_bh')
pval = -np.log10(pval)
corr_pvals = -np.log10(corr_pvals)

if '--bothsided' in sys.argv:
    pval = np.sign(tt)*pval
    corr_pvals = np.sign(tt)*corr_pvals


np.savetxt(os.path.splitext(sys.argv[1])[0]+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+'_tstats.txt', np.concatenate([corr[:,[0]], np.around(tt,3).reshape(-1,1), np.around(pval,3).reshape(-1,1), np.around(corr_pvals,3).reshape(-1,1)], axis = 1).astype(str), fmt = '%s', header = 'Gene Tstat log10pvalue BH_corrected')



