# plots fraction of drivertypes per gene against correlation to observed data
#python3 ../plot_drivertype.py ALLgenes_ism_attributions_driversfw_types.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list -weighted ALLgenes_ism_attributions_driversfw.txt -3

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from functools import reduce
import seaborn as sns
from scipy.stats import mannwhitneyu as mwu


# Read combined drivertype file that contains drivers for all genes
counts = np.genfromtxt(sys.argv[1], dtype = str)
allctypes = counts[:,-2:].astype(float)
bitypes = np.zeros((len(counts),2))
# Assign is-supported and is unsupported driver
for c, ctypes in enumerate(allctypes):
    if (ctypes[0]<0 and ctypes[1]<0) or (ctypes[0]>0 and ctypes[1]>0):
        bitypes[c,0] = 1
    elif ctypes[0]>0 and ctypes[1]<0:
        bitypes[c,1] = 1
    elif ctypes[0]<0 and ctypes[1]>0:
        bitypes[c,1] = 1
# Read prediction values
predictions = np.genfromtxt(sys.argv[2], dtype = str)
# Only investigate genes that can be approximated with the sum of ISMs > 0.2 to prediction 
approx_control = np.genfromtxt(sys.argv[3], dtype = str)

# Sort files to align
commons = np.sort(reduce(np.intersect1d, [counts[:,0], predictions[:,0], approx_control[:,0]]))
a_, sort = np.unique(approx_control[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
approx_control = approx_control[sort]
a_, sort = np.unique(predictions[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
predictions = predictions[sort]

# Select subset of genes
if '--genelist' in sys.argv:
    genelist = np.genfromtxt(sys.argv[sys.argv.index('--genelist')+1], dtype = str)
    print(np.array_equal(approx_control[:,0], commons))
    mask = np.isin(commons, genelist)
    approx_control, predictions, commons = approx_control[mask], predictions[mask], commons[mask]

# Weight drivers by the attribution to the prediction for mean within each gene
weights = np.ones(len(counts))
if '--weighted' in sys.argv:
    wfile = np.genfromtxt(sys.argv[sys.argv.index('--weighted')+1], dtype = str)
    wl = int(sys.argv[sys.argv.index('--weighted')+2])
    cnames = np.array([cnt[0]+'_'+cnt[1] for cnt in counts])
    wnames = np.array([cnt[0]+'_'+cnt[1] for cnt in wfile])
    csort = np.argsort(cnames)[np.isin(np.sort(cnames), wnames)]
    counts, types, bitypes = counts[csort], types[csort], bitypes[csort]
    weights = wfile[np.argsort(wnames)[np.isin(np.sort(wnames), cnames)],wl].astype(float)

# Generate mean per gene from individual drivers
bipgcounts = np.zeros((len(commons),2))
keep = np.zeros(len(commons)) == 0
for c, com in enumerate(commons):
    if com in counts[:,0]:
        cmask = counts[:,0]==com
        bipgcounts[c] = np.sum(weights[cmask][:,None]*bitypes[cmask] , axis = 0)/np.sum(weights[cmask])
    else:
        keep[c] = False

approx_control, predictions, commons, bipgcounts = approx_control[keep], predictions[keep], commons[keep], bipgcounts[keep]

#Only investigate genes that can be approximated with the sum of ISMs > 0.2 to prediction  
control = approx_control[:,1].astype(float) > 0.2
approx_control, predictions, commons, pgcounts, bipgcounts = approx_control[control], predictions[control], commons[control], pgcounts[control], bipgcounts[control]

print('Commons fine', np.array_equal(commons, approx_control[:,0]))

# Divide sets into sets with significant positive and negative correlation 
bpos = predictions[:,1].astype(float) > 0.2
bneg = predictions[:,1].astype(float) < -0.2

# Assign p-values to differences between two sets and replace value with * symbols
pvalcut = np.array([0.05, 0.01,0.001])
pvalsign = np.array(['*', '**', '***'])

# Print some stats
print('Positive')
for bpo in np.where(bpos)[0]:
    print(commons[bpo], bipgcounts[bpo])

print('\nNegative')
for bpo in np.where(bneg)[0]:
    print(commons[bpo], bipgcounts[bpo])

# Generate figures for fraction of supported and unsupported drivers
muttype = ['supported', 'unsupported']
for t, mut in enumerate(muttype):
    fig = plt.figure(figsize = (1.5,4))
    ax = fig.add_subplot(111)
    parts = ax.violinplot([bipgcounts[bpos,t], bipgcounts[bneg,t]], widths = 0.95, showmeans=False, showmedians=False, showextrema=True)
    fcolors = ['forestgreen', 'indigo']
    for p, pc in enumerate(parts['bodies']):
        pc.set_facecolor(fcolors[p])
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    for partname in ('cbars','cmins','cmaxes'):
        vp = parts[partname]
        vp.set_edgecolor('k')
    quartile1, medians, quartile3 = [np.percentile(bipgcounts[bpos,t], [25]),np.percentile(bipgcounts[bneg,t], [25])], [np.percentile(bipgcounts[bpos,t], [50]),np.percentile(bipgcounts[bneg,t], [50])],[np.percentile(bipgcounts[bpos,t], [75]),np.percentile(bipgcounts[bneg,t], [75])]
    ax.scatter([1,2], medians, marker='o', color='grey', s=30, zorder=3)
    ax.vlines([1,2], quartile1, quartile3, color='k', linestyle='-', lw=10)
    pval = mwu(bipgcounts[bpos,t], bipgcounts[bneg,t])[1]*len(muttype)/2
    ps = np.where(pvalcut>pval)[0]
    if len(ps) > 0:
        ax.plot([1,2], [1.1, 1.1], c = 'k')
        ax.plot([1,1], [1.07,1.1], c = 'k')
        ax.plot([2,2], [1.07, 1.1], c = 'k')
        ax.text(1.5, 1.11, pvalsign[ps[-1]], ha = 'center', va = 'bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([1,2])
    ax.set_xticklabels(['Positive', 'Negative'], rotation = 45)
    ax.set_xlabel('Correlation to Enformer')
    ax.set_ylabel('Fraction of '+mut+' drivers')
    fig.savefig(os.path.splitext(sys.argv[1])[0]+'_boxplot_negvspos'+mut+'.svg', transparent = True, dpi = 450, bbox_inches = 'tight')
    
    
    figx = plt.figure(figsize = (3.,3.5))
    axx = figx.add_subplot(111)
    axx.spines['top'].set_visible(False)
    axx.spines['right'].set_visible(False)
    axx.plot([0,0],[0,1], color = 'grey', zorder = -1)
    axx.scatter(predictions[bpos,1].astype(float), bipgcounts[bpos,t], color = 'forestgreen')
    axx.scatter(predictions[bneg,1].astype(float), bipgcounts[bneg,t], color = 'indigo')
    axx.scatter(predictions[~bneg*~bpos,1].astype(float), bipgcounts[~bneg*~bpos,t], color = 'grey')
    axx.set_xlabel('Correlation to Enformer')
    axx.set_ylabel('Fraction of '+mut+' drivers')

    figx.savefig(os.path.splitext(sys.argv[1])[0]+'_scatter_negvspos'+mut+'.jpg', transparent = True, dpi = 450, bbox_inches = 'tight')


