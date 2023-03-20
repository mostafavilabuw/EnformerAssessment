# Performs statistical test for sequence patterns to be enriched for unsupported or supported snp
#python3 ../plot_common_motifs.py ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.25 --savefig ALLgenes_ism_attributions_driversfw_refseq_winsize13_clust_ms5-1_complete0.2

import numpy as np
import sys, os
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.multitest import multipletests
import logomaker as lm

def read_meme(mfile):
    obj = open(mfile, 'r').readlines()
    pwms = []
    pwmnames = []
    pwm = []
    t = -1000
    for l, line in enumerate(obj):
        line = line.strip().split()
        if len(line) > 0:
            if line[0] == 'ALPHABET=':
                nts = list(line[1])
            if line[0] == 'MOTIF':
                pwmnames.append(line[1].split('_')[1])
                t = 0
            t += 1
            if t > 2:
                pwm.append(line)
        else:
            if len(pwm) > 0:
                pwm = np.array(pwm, dtype = float)
                #pwm = np.log2(pwm/0.25)
                #pwm[pwm<0] = 0
                pwms.append(pd.DataFrame(data = np.array(pwm, dtype = float), columns = nts))
                pwm = []
    if len(pwm) > 0:
        pwm = np.array(pwm, dtype = float)
        #pwm = np.log2(pwm/0.25)
        #pwm[pwm<0] = 0
        pwms.append(pd.DataFrame(data = pwm, columns = nts))
    return pwms, pwmnames
        
# Read driver type
drivertype = np.genfromtxt(sys.argv[1], dtype = str) # to select and assign type

# select gene set
if '--geneset' in sys.argv:
    geneset = np.genfromtxt(sys.argv[sys.argv.index('--geneset')+1], dtype = str)
    dmask = np.isin(drivertype[:,0], geneset)
    drivertype = drivertype[dmask]

allctypes = drivertype[:,-2:].astype(float)
drivers = np.array([d[0]+d[1] for d in drivertype])
mask = np.argsort(drivers)
drivers, allctypes = drivers[mask], allctypes[mask]
bitypes = np.zeros((len(drivers),2))
for c, ctypes in enumerate(allctypes):
    #print(ctypes)
    if (ctypes[0]<0 and ctypes[1]<0) or (ctypes[0]>0 and ctypes[1]>0):
        bitypes[c,0] = 1
    elif ctypes[0]>0 and ctypes[1]<0:
        bitypes[c,1] = 1
    elif ctypes[0]<0 and ctypes[1]>0:
        bitypes[c,1] = 1

# Start of file names that contain the information about the clustered sequences
clufile = sys.argv[2]

if '--savefig' in sys.argv:
    outname = sys.argv[sys.argv.index('--savefig')+1]

# Read pwms of clusters
pwms, pwmnames = read_meme(clufile+'_clusterpwms.txt')

# read cluster assignments
clusters = np.genfromtxt(clufile+'_clusters.txt', dtype = str)
clusternames, clusters = clusters[:,:2], clusters[:,2]
clusternames = np.array([d[0]+d[1] for d in clusternames])

# read location of snps in these pwms
snploc = np.genfromtxt(clufile+'_locclpwms.txt', dtype = str)
snpname, loc = snploc[:,:2], snploc[:,2].astype(int)
snpname = np.array([d[0]+d[1] for d in snpname])

if not np.array_equal(snpname, clusternames):
    print('CLUSTERFiles are different')
    sys.exit()
# align clusternames and snpnames to drivers
mask = np.argsort(snpname)[np.isin(np.sort(snpname), drivers)]
clusters, clusternames, loc, snpname = clusters[mask], clusternames[mask], loc[mask], snpname[mask]

#align drivers to clusternames
mask = np.argsort(drivers)[np.isin(np.sort(drivers), snpname)]
drivers, types, bitypes = drivers[mask], types[mask], bitypes[mask]
print(np.array_equal(snpname, drivers), len(snpname), len(drivers))

cnames, cnum = np.unique(clusters, return_counts=True)
print('N clusters', len(cnames), 'in', len(clusters))
from scipy.stats import fisher_exact
pvalsneg = []
pvalspos = []
totalsize = []
# for each cluster, use fisher exact test to compute enrichment of supported and unsupported snps
for c, cna in enumerate(cnames):
    table = [[np.sum((clusters == cna) * (bitypes[:,1] == 1)), np.sum((clusters == cna) * (bitypes[:,1] == 0))],[np.sum((clusters != cna) * (bitypes[:,1] == 1)), np.sum((clusters != cna) * (bitypes[:,1] == 0))]]
    oddspos, p_valuepos = fisher_exact(table, 'less') # pvalue that positives are enriched
    oddsneg, p_valueneg = fisher_exact(table, 'greater') # pvalue that negatives are enriched
    pvalsneg.append(p_valueneg)
    pvalspos.append(p_valuepos)
    totalsize.append(np.sum(clusters == cna))


totalsize = np.array(totalsize)
pvalsneg = np.array(pvalsneg)
pvalspos = np.array(pvalspos)

# adjsut p-value with benjamini hochberg
ylabel = 'p_value'
if '--adjust_pvalue' in sys.argv:
    is_, pvalsneg, a_, b_ = multipletests(pvalsneg, method='fdr_bh')
    is_, pvalspos, a_, b_ = multipletests(pvalspos, method='fdr_bh')
    ylabel = 'p_value_corrected_BH'
pvalsneg = -np.log10(pvalsneg)
pvalspos = -np.log10(pvalspos)

sortneg = np.argsort(pvalsneg)
sortpos = np.argsort(-pvalspos)

fig = plt.figure(figsize = (6,4))
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.scatter(np.arange(len(pvalsneg)), -pvalsneg[sortneg], s = totalsize[sortneg] * 5, label = 'unsupported')
ax.scatter(np.arange(len(pvalspos)), pvalspos[sortpos], s = totalsize[sortpos] * 5, label = 'supported')
ax.set_ylabel(ylabel)
ax.set_xlabel('Motif cluster')
ax.plot([0,len(pvalsneg)], [0,0],color = 'k')
ax.plot([0,len(pvalsneg)], [-np.log10(0.01),-np.log10(0.01)],color = 'grey',ls = '--')
ax.plot([0,len(pvalsneg)], [-np.log10(0.05),-np.log10(0.05)],color = 'grey',ls = '--')
ax.plot([0,len(pvalsneg)], [-np.log10(0.05/len(pvalsneg)),-np.log10(0.05/len(pvalsneg))],color = 'red',ls = '--')
ax.plot([0,len(pvalsneg)], [np.log10(0.05),np.log10(0.05)],color = 'grey',ls = '--')
ax.plot([0,len(pvalsneg)], [np.log10(0.01),np.log10(0.01)],color = 'grey',ls = '--')
ax.plot([0,len(pvalsneg)], [np.log10(0.05/len(pvalsneg)),np.log10(0.05/len(pvalsneg))],color = 'red',ls = '--')
ax.set_yticks([-np.log10(0.05/len(pvalsneg)), -np.log10(0.01),-np.log10(0.05),0,np.log10(0.05),np.log10(0.01),np.log10(0.05/len(pvalsneg))])
ax.set_yticklabels(['FDR 0.05', '0.01', '0.05', '0', '0.05', '0.01', 'FDR 0.05'])
ax.legend()
if '--savefig' in sys.argv:
    fig.savefig(outname + '_clustenrich.jpg', dpi = 200, bbox_inches = 'tight')
else:
    plt.show()

fig2 = plt.figure(figsize = (4,4))
# plot distribution of clustersizes
ax2 = fig2.add_subplot(111)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
bins = np.arange(np.amax(cnum)+2)
ax2.hist(cnum, bins = bins)
ax2.set_xlabel('Cluster size')
if '--savefig' in sys.argv:
    fig2.savefig(outname + '_clustsizedist.jpg', dpi = 200, bbox_inches = 'tight')
else:
    plt.show()


# plot pwms of most significant clusters
for s in sortneg[::-1][:20]:
    clu = cnames[s]
    clumask = clusters == clu
    pwm = pwms[pwmnames.index(clu)]
    ploc = loc[clumask]
    snpdir = bitypes[clumask,0] - bitypes[clumask, 1]
    dcounts = []
    for sd in [-1,1]:
        dcounts.append(np.unique(ploc[snpdir == sd], return_counts = True))
    figp = plt.figure(figsize = (len(pwm)*0.5,2))
    axp = figp.add_subplot(111)
    axp.spines['top'].set_visible(False)
    axp.spines['right'].set_visible(False)
    lm.Logo(pwm, ax = axp)
    scale = -0.05
    axp.set_title('logp; '+str(round(pvalsneg[s],3))+'; numvar:'+str(int(cnum[s])))
    axp.bar(dcounts[0][0],scale* dcounts[0][1], label = 'False', width = 0.4)
    bottom = np.zeros(len(dcounts[1][0]))
    bottom[np.isin(dcounts[1][0], dcounts[0][0])] = dcounts[0][1][np.isin(dcounts[0][0],dcounts[1][0])]
    axp.bar(dcounts[1][0], scale*dcounts[1][1], bottom = scale*bottom, label = 'Correct', width = 0.4)
    min0, min1 = 0,0
    if len(dcounts[0][1]) > 0:
        min0 = np.amin(scale*dcounts[0][1])
    if len(dcounts[1][1]) > 0:
        min1 = np.amin(scale*dcounts[1][1]+bottom)

    axp.set_ylim([min(min0,min1)+scale, 1])
    axp.set_yticks([min(min0,min1), 0.5,1])
    axp.set_yticklabels([int(min(min0,min1)/scale), 0.5, 1])
    if '--savefig' in sys.argv:
        figp.savefig(outname + '_clust'+str(clu)+'.jpg', dpi = 200, bbox_inches = 'tight')
        print(outname + '_clust'+str(clu)+'.jpg')
    else:
        plt.show()
    plt.close()

