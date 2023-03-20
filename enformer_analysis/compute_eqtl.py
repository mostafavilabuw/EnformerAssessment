# From observed gene expression across individuals and thei genoptypes compute the eqtl value and correlation with observed gene expression

import numpy as np
import sys, os
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

# Set of genes for which eqtls will be computed
tgenes = np.genfromtxt(sys.argv[1], dtype = str)

# load observed gene exppression
obsfile = open('../Observed_gene_expressionENSG.tsv','r').readlines()
inds = []
obs_exp = []
for l, line in enumerate(obsfile):
    line = line.strip().split()
    if l == 0:
        genes = line
    else:
        inds.append(line[0])
        obs_exp.append(line[1:])

genes = np.array(genes)
inds = np.array(inds)
obs_exp = np.array(obs_exp, dtype = float)

for i, gene in enumerate(tgenes):
    print(gene)
    g = list(genes).index(gene)
    obs_e = obs_exp[:,g]
    
    # Load SNV file 
    snpfile = np.load('../variant_info_100k/'+gene+'snp_info.npz' )
    snp_info = snpfile['snps'].astype(float)
    snp_name = snpfile['rows'].astype(str)
    ind_names = snpfile['columns'].astype(str)
    
    sortobs = np.argsort(inds)[np.isin(np.sort(inds), ind_names)]
    obs_e = obs_e[sortobs]

    sortvar = np.argsort(ind_names)[np.isin(np.sort(ind_names), inds)]
    snp_info = snp_info[:, sortvar]
    ind_names = ind_names[sortvar]
    # check sorting of individuals
    print(np.array_equal(inds[sortobs], ind_names))

    obj = open(gene+'_eqtl.txt', 'w')
    obj2 = open(gene+'_corr.txt', 'w')
    
    lr = LinearRegression(fit_intercept=True)
    corr = []
    eqtl = []
    for s, sn in enumerate(snp_name):
        pear = pearsonr(snp_info[s], obs_e)[0]
        lr = lr.fit(snp_info[s].reshape(-1,1),obs_e)
        coef_ = lr.coef_
        obj2.write(sn+' '+str(round(pear,4))+'\n')
        obj.write(sn+' '+str(round(coef_[0],6))+'\n')
        corr.append(pear)
        eqtl.append(coef_[0])
    argmax = np.argmax(np.absolute(np.nan_to_num(corr)))
    print(corr[argmax], eqtl[argmax])
    argmax = np.argmax(np.absolute(np.nan_to_num(eqtl)))
    print(corr[argmax], eqtl[argmax])
    obj.close()
    obj2.close()
