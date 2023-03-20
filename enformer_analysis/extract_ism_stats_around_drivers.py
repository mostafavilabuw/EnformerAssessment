# Takes ISM attribution windows around the drivers and compares to "global" attributions stats from the entire sequence
# Determines standardized effect of base change
# Determine size of susequent motifs at drivers with certain impact
# Determine distance to base with certain impact

import numpy as np
import matplotlib.pyplot as plt
import sys, os

# Read drivers
drivers = np.genfromtxt(sys.argv[1], dtype = str) # file with all drivers

# Select subset of genes
if '--geneset' in sys.argv:
    geneset = np.genfromtxt(sys.argv[sys.argv.index('--geneset')+1], dtype = str)
    drivers = drivers[np.isin(drivers[:,0], geneset)]

genes = drivers[:,0]
gloc = drivers[:,1].astype(int)

# Read in global stats for isms
#globalstds = np.genfromtxt('baseline_ism/Baseline_stats.txt', dtype = str)
globalstds = np.genfromtxt('../TSS_ISM/tss_attribution_stats.txt', dtype = str)

nts = list('ACGT')
ntsar = np.array(nts)

ar = []
for g, gene in enumerate(genes):
    # load the attributions around the driver snps at the reference sequence
    std = float(globalstds[list(globalstds[:,0]).index(gene),2])
    attatreffile = np.load('from_ref/'+gene+'_'+str(gloc[g])+'_attribution.npz')
    attatref, seqatref = attatreffile['attribution'], str(attatreffile['seq'])
    attatref/=std
    attatref = -np.sum(attatref, axis = 1)/3
    # load the attributions around the driver snps at the variant sequence
    attatvarfile = np.load('from_main_var/'+gene+'_'+str(gloc[g])+'_attribution.npz')
    attatvar, seqatvar = attatvarfile['attribution'], list(str(attatreffile['seq']))
    attatvar/=std
    attatvar = -np.sum(attatvar, axis = 1)/3

    # check check direction to each other of zscore and motifs
    if '--checkmax' in sys.argv:
        thresh = np.array([0.05,0.1,0.2,0.5])
    else:
        thresh = np.array([1.64,1.96,2.58,3.29])
    print(thresh)
    loc = int(len(attatref)/2)
    latt = len(attatref)
    
    score = attatref[loc]
    size = (thresh<=abs(score)).astype(int)
    dist = (thresh>abs(score)).astype(int)
    checksizeneg = thresh <= abs(score)
    checksizepos = thresh <= abs(score)
    checkdist = thresh > abs(score)
    i = 1
    # check the size of subsequent bases with thresh or how far first base is from snp with thresh
    while True:
        #print(loc-i,loc+j, score, checkdist, checksizeneg, checksizepos, attribution[loc+i], attribution[loc-i], dist)
        if loc + i < latt:
            checksizepos = checksizepos * (thresh <= (np.sign(score) * attatref[loc+i]))
            size += checksizepos.astype(int)
            checkdist = checkdist * (thresh > abs(attatref[loc+i]))
        if loc - i >=0:
            checksizeneg = checksizeneg * (thresh <= (np.sign(score) * attatref[loc-i]))
            size += checksizeneg.astype(int)
            checkdist = checkdist * (thresh > abs(attatref[loc-i]))
        if loc - i <=0 or loc + i >= latt-1:
            break
        dist += checkdist.astype(int)
        if not checksizepos.any() and not checksizeneg.any() and not checkdist.any():
            break
        i += 1
        
    distscore = np.argmax(np.absolute(np.array([attatref[loc+dist], attatref[loc-dist]])), axis = 0)
    distscore = np.array([attatref[loc+dist], attatref[loc-dist]])[distscore,np.arange(4,dtype = int)]
    add = np.concatenate([[gene, str(gloc[g]), str(score)], size.astype(str), dist.astype(str), distscore.astype(str)])

    score = attatvar[loc]
    size = (thresh<=abs(score)).astype(int)
    dist = (thresh>abs(score)).astype(int)
    checksizeneg = thresh <= abs(score)
    checksizepos = thresh <= abs(score)
    checkdist = thresh > abs(score)
    i = 1
    while True:
        #print(loc-i,loc+j, score, checkdist, checksizeneg, checksizepos, attribution[loc+i], attribution[loc-i], dist)
        if loc + i < latt:
            checksizepos = checksizepos * (thresh <= (np.sign(score) * attatvar[loc+i]))
            size += checksizepos.astype(int)
            checkdist = checkdist * (thresh > abs(attatvar[loc+i]))
        if loc - i >=0:
            checksizeneg = checksizeneg * (thresh <= (np.sign(score) * attatvar[loc-i]))
            size += checksizeneg.astype(int)
            checkdist = checkdist * (thresh > abs(attatvar[loc-i]))
        if loc - i <=0 or loc + i >= latt-1:
            break
        dist += checkdist.astype(int)
        if not checksizepos.any() and not checksizeneg.any() and not checkdist.any():
            break
        i += 1

    distscore = np.argmax(np.absolute(np.array([attatvar[loc+dist], attatvar[loc-dist]])), axis = 0)
    distscore = np.array([attatvar[loc+dist], attatvar[loc-dist]])[distscore,np.arange(4,dtype = int)]
    add = np.append(add,np.concatenate([[str(score)], size.astype(str), dist.astype(str), distscore.astype(str)]))

    ar.append(add)
    print(ar[-1])

ar = np.array(ar)
print(os.path.splitext(os.path.split(sys.argv[1])[1])[0] + '_ism_significance_stats.txt')
np.savetxt(os.path.splitext(os.path.split(sys.argv[1])[1])[0] + '_ism_significance_stats.txt', ar.astype(str), fmt='%s')


