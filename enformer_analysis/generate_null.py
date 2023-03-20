# Generates random Null distribution for absolute correlation between observed and predicted gene expression from set of SNVs for each individual
# For each gene, SNVs get random attribution assigned and the sum of attributions times the SNVs in a person's genome is the predicted gene expression.
# This random process is performed several times and the mean and std of the absolute correlation are computed as a random Null.
# Basic usage:
# python3 generate_null.py


import numpy as np
import sys, os
import glob
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def correlation(y1, y2, axis = 1, ctype = 'pearson', distance = False):
    if ctype == 'spearman':
        y1, y2 = np.argsort(np.argsort(y1, axis = axis), axis = axis), np.argsort(np.argsort(y2,axis = axis), axis = axis)
    if ctype != 'cosine':
        mean1, mean2 = np.mean(y1, axis = axis), np.mean(y2, axis = axis)
        y1mean, y2mean = y1-np.expand_dims(mean1,axis = axis), y2-np.expand_dims(mean2,axis = axis)
    else:
        y1mean, y2mean = y1, y2
    n1, n2 = np.sqrt(np.sum(y1mean*y1mean, axis = axis)), np.sqrt(np.sum(y2mean*y2mean, axis = axis))
    n12 = n1*n2
    y12 = np.sum(y1mean*y2mean, axis = axis)
    if isinstance(y12, float):
        if n12/max(n1,n2) < 1e-8:
            n12, y12 = 1., -1.
        else:
            y12[n12/np.amax(np.array([n1,n2]),axis = 0) < 1e-8] = -1.
            n12[n12/np.amax(np.array([n1,n2]),axis = 0) < 1e-8] = 1
    corout = y12/n12
    if distance:
        corout = 1.-corout
    return np.around(corout,4)


# Find all snp_info.npz files
sfiles = np.sort(glob.glob('ENSG*snp_info.npz'))
# Load ID of indivdiuals
inds = np.load(sfiles[0])['columns']
# Make list of gene names
genes = np.array([s.split('snp_info')[0] for s in sfiles])

# Load Observed gene expression data
obs = np.genfromtxt('Observed_gene_expressionENSG.tsv', dtype = str, skip_header = 1)
obsgenes = np.array(open('Observed_gene_expressionENSG.tsv', 'r').readline().strip().split('\t'))
obsind, obs = obs[:,0], obs[:,1:].astype(float)
# Sort observed expression to given genes
osort = np.argsort(obsgenes)[np.isin(np.sort(obsgenes), genes)]
obsgenes, obs = obsgenes[osort], obs[:,osort]
# sort observed expression to given individuals
osort = np.argsort(obsind)[np.isin(np.sort(obsind), inds)]
obsind, obs = obsind[osort], obs[osort]
# sort given inds to to observed expression inds
isort = np.argsort(inds)[np.isin(np.sort(inds), obsind)]
inds = inds[isort]
# check if sorting was successful 
print(np.array_equal(obsind, inds), len(obsind), len(inds))

np.random.seed(1)
nstat = 50

# Second file that compares expression to predicted values can be generated
if '--savemeanexp' in sys.argv:
    enf = np.genfromtxt('Enformer_predictionsENSG.tsv', dtype = str, skip_header = 1)
    enfgenes = np.array(open('Enformer_predictionsENSG.tsv', 'r').readline().strip().split('\t'))
    enfind, enf = enf[:,0], enf[:,1:].astype(float)
    esort = np.argsort(enfgenes)[np.isin(np.sort(enfgenes), genes)]
    enfgenes, enf = enfgenes[esort], enf[:,esort]
    esort = np.argsort(enfind)[np.isin(np.sort(enfind), inds)]
    enfind, enf = enfind[esort], enf[esort]    
    print('enformer', np.array_equal(enfind, inds))


corrs = [] # mean std n
exps = []
for g, gene in enumerate(genes):
    
    snpfile = np.load(gene+'snp_info.npz')
    indv = snpfile['columns'][isort]
    if np.array_equal(indv, obsind):
        snp, snames = snpfile['snps'], snpfile['rows']
        snp = snp[:,isort]

        rat= np.random.normal(size = (len(snames), nstat))
        # Different option for Null with all SNV weights being positive
        if '--positive' in sys.argv:
            rat = np.absolute(rat)
        ex = np.sum(snp.T[...,None]*rat[None,...], axis = 1)
        pears = correlation(ex, obs[:,np.where(obsgenes == gene)[0][[0]]], axis = 0)
        pears = np.absolute(np.nan_to_num(pears))
        corrs.append([gene, np.around(np.mean(pears),3), np.around(np.std(pears),5), nstat])
        print(corrs[-1])
        if '--savemeanexp' in sys.argv:
            pears = correlation(ex, enf[:,np.where(enfgenes == gene)[0][[0]]], axis = 0)
            pears = np.absolute(np.nan_to_num(pears))
            exps.append([gene, np.around(np.mean(pears),3), np.around(np.std(pears),5), nstat])
            print(exps[-1])

corrs = np.array(corrs, dtype = str)
if '--savemeanexp' in sys.argv:
    exps = np.array(exps, dtype  = str)

if '--positive' in sys.argv:
    np.savetxt('GeneSpecific_CorrelationtoObsRandomPosNull.txt', corrs, header = 'Gene MeanR StdR N', fmt = '%s')
    if '--savemeanexp' in sys.argv:
        np.savetxt('GeneSpecific_CorrelationtoObsRandomPosNulltoEnf.txt', exps, header = 'Gene MeanR StdR N', fmt = '%s')
else:
    np.savetxt('GeneSpecific_CorrelationtoObsRandomNull.txt', corrs, header = 'Gene MeanR StdR N', fmt = '%s')
    if '--savemeanexp' in sys.argv:
        np.savetxt('GeneSpecific_CorrelationtoObsRandomNulltoEnf.txt', exps, header = 'Gene MeanR StdR N', fmt = '%s')



