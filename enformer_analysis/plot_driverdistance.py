# Plot histogram of the distances between drivers and the TSS
# 
#python3 ../plot_driverdistance.py DistancetoTSS_ism_attributions_driversbw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --plot_main_in_all

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from functools import reduce

# Read Distances to TSS
counts = np.genfromtxt(sys.argv[1], dtype = str)
# Read correletions of predictions
predictions = np.genfromtxt(sys.argv[2], dtype = str)
# Read in correlation of sum of attributions to prediction with full model
approx_control = np.genfromtxt(sys.argv[3], dtype = str)

# sort files to match
commons = np.sort(reduce(np.intersect1d, [counts[:,0], predictions[:,0], approx_control[:,0]]))
a_, sort = np.unique(approx_control[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
approx_control = approx_control[sort]
a_, sort = np.unique(predictions[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
predictions = predictions[sort]
counts = counts[np.argsort(counts[:,0])]

# check if sorting worked
print(np.array_equal(predictions[:,0], approx_control[:,0]), np.array_equal(commons, predictions[:,0]), len(predictions))

# Select subset of genes
if '--genelist' in sys.argv:
    genelist = np.genfromtxt(sys.argv[sys.argv.index('--genelist')+1], dtype =str)
    mask = np.isin(predictions[:,0], genelist)
    predictions, approx_control, commons = predictions[mask], approx_control[mask], commons[mask]

# Only work with genes that can be linearly approximated
control = approx_control[:,1].astype(float) > 0.2
approx_control, predictions, commons = approx_control[control], predictions[control], commons[control]
print('Commons fine', np.array_equal(commons, approx_control[:,0]))

counts = counts[np.isin(counts[:,0], approx_control[:,0])]
bcountmask = []
for c, com in enumerate(commons):
    isgene = np.where(counts[:,0] == com)[0]
    if len(isgene) > 0:
        bcountmask.append(isgene[np.argmax(counts[isgene,-1].astype(float))])
bcounts = counts[bcountmask]

# Split genes in genes with postive and negative correlation to observed data
bpos = predictions[:,1].astype(float) > 0.2
bneg = predictions[:,1].astype(float) < -0.2
pos = np.isin(counts[:,0], commons[bpos])
neg = np.isin(counts[:,0], commons[bneg])


# Generate histograms
fig = plt.figure(figsize = (6,4))
axp = fig.add_subplot(121)
axn = fig.add_subplot(122)

bins = np.concatenate([-np.arange(10000,10000,10000)[::-1],-np.arange(2000,10000,2000)[::-1],-np.arange(500,2000,250)[::-1],-np.arange(0,500,100)[::-1], np.arange(0,500,100),np.arange(500,2000,250),np.arange(2000,10000,2000), np.arange(10000,10000,10000)])

bins = np.arange(-70000,70000,1000)

a,b,c = axp.hist(counts[pos,2].astype(float), color = 'forestgreen', alpha = 0.9, bins = bins)
axp.set_title('Positively correlated genes')
axp.spines['top'].set_visible(False)
axp.spines['right'].set_visible(False)
axp.set_xlabel('Distance to TSS')
axp.set_ylabel('Number drivers')
if '--plot_main_in_all' in sys.argv:
    axp.hist(bcounts[bpos,2].astype(float), color = 'orange', bins = bins, alpha = 0.9, label = 'Main\ndrivers')
    axp.legend(loc = 'upper right', fontsize = 'small')
axp.plot([0,0],[0,np.amax(a)], c = 'goldenrod', ls = '--')
#axp.set_xscale('symlog')


a,b,c = axn.hist(counts[neg,2].astype(float), color = 'indigo', alpha = 0.9, bins = bins)
axn.set_title('Negatively correlated genes')
axn.spines['top'].set_visible(False)
axn.spines['right'].set_visible(False)
axn.set_xlabel('Distance to TSS')
axn.set_ylabel('Number drivers')

# Show main drivers in different color
if '--plot_main_in_all' in sys.argv:
    axn.hist(bcounts[bneg,2].astype(float), color = 'mediumvioletred', bins = bins, alpha = 0.9, label = 'Main\ndrivers')
    axn.legend(loc = 'upper right', fontsize = 'small')
axn.plot([0,0],[0,np.amax(a)], c = 'goldenrod', ls = '--')
#axn.set_xscale('symlog')


# Generate extra figure main drivers
figc = plt.figure(figsize = (6,4))
axpc = figc.add_subplot(121)
axnc = figc.add_subplot(122)

axpc.hist(bcounts[bpos,2].astype(float), color = 'forestgreen', bins = bins, alpha = 0.5)
axpc.set_title('Positively correlated genes')
axpc.spines['top'].set_visible(False)
axpc.spines['right'].set_visible(False)
axpc.set_xlabel('Distance to TSS')
axpc.set_ylabel('Number main drivers')

axnc.hist(bcounts[bneg,2].astype(float), color = 'indigo', bins = bins, alpha = 0.5)
axnc.set_title('Negatively correlated genes')
axnc.spines['top'].set_visible(False)
axnc.spines['right'].set_visible(False)
axnc.set_xlabel('Distance to TSS')
axnc.set_ylabel('Number main drivers')

fig.savefig(os.path.splitext(sys.argv[1])[0]+'_distribution.jpg', dpi = 450, bbox_inches = 'tight')
axp.set_xlim([-19999,19999])
axn.set_xlim([-19999,19999])
fig.savefig(os.path.splitext(sys.argv[1])[0]+'_distributionzoom.jpg', dpi = 450, bbox_inches = 'tight')
figc.savefig(os.path.splitext(sys.argv[1])[0]+'_maindistribution.jpg', dpi = 450, bbox_inches = 'tight')

