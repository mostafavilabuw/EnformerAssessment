# plot the distribution of the number of drivers per investigated genes
#python3 ../plot_drivercounts.py Counts_ism_attributions_driversfw.txt ../Prediction_correlationsCageAdultBrain.txt ALL_genes_ism_attributions_sumpersonal_mp_vs_Enformer_predictions_CageAdultBrain.txt --genelist ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list --split_sets

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from functools import reduce

# Read driver counts
counts = np.genfromtxt(sys.argv[1], dtype = str)
# Read correlation of predictions to observed
predictions = np.genfromtxt(sys.argv[2], dtype = str)
# Read correlation between prediction and linear approximation
approx_control = np.genfromtxt(sys.argv[3], dtype = str)

# sort files to match rows
commons = reduce(np.intersect1d, [counts[:,0], predictions[:,0], approx_control[:,0]])
a_, sort = np.unique(approx_control[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
approx_control = approx_control[sort]
a_, sort = np.unique(predictions[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
predictions = predictions[sort]
a_, sort = np.unique(counts[:,0], return_index = True)
sort = sort[np.isin(a_, commons)]
counts = counts[sort]

# Check if sorting was successful
if not np.array_equal(counts[:,0], approx_control[:,0]) or not np.array_equal(counts[:,0], predictions[:,0]): 
    print('Files not alinged')
    print(np.array_equal(counts[:,0], approx_control[:,0]), np.array_equal(counts[:,0], predictions[:,0]), len(predictions))
    sys.exit()
    
# select subset of genes
if '--genelist' in sys.argv:
    genelist = np.genfromtxt(sys.argv[sys.argv.index('--genelist')+1], dtype =str)
    mask = np.isin(counts[:,0], genelist)
    counts, predictions, approx_control = counts[mask], predictions[mask], approx_control[mask]

control = approx_control[:,1].astype(float) > 0.2
counts, predictions = counts[control], predictions[control]


# split into positively and negatively correlated genes
if '--split_sets' in sys.argv:
    pos = predictions[:,1].astype(float) > 0.2
    neg = predictions[:,1].astype(float) < -0.2
else:
    pos = np.ones(len(predictions)) ==1
    neg = np.ones(len(predictions)) ==1



fig = plt.figure(figsize = (8,4))
axp = fig.add_subplot(121)
bins = np.arange(np.amax(counts[pos|neg,1].astype(float))+1)-0.5
a,b,c = axp.hist(counts[pos,1].astype(float), color = 'navy', alpha = 0.5, bins = bins)
print(a)
print(b)
if '--split_sets' in sys.argv:
    axp.set_title('Positive correlated')
axp.spines['top'].set_visible(False)
axp.spines['right'].set_visible(False)
axp.set_xlabel('Number drivers')
axp.set_ylabel('Number genes')
axp.set_xlim([0,20])

if '--split_sets' in sys.argv:
    axn = fig.add_subplot(122)
    a,b,c = axn.hist(counts[neg,1].astype(float), color = 'navy', alpha = 0.5, bins = bins)
    axn.set_title('Negative correlated')
    axn.spines['top'].set_visible(False)
    axn.spines['right'].set_visible(False)
    axn.set_xlabel('Number drivers')
    axn.set_ylabel('Number genes')
    axn.set_xlim([0,20])


figc = plt.figure(figsize = (8,4))
axpc = figc.add_subplot(121)

a_, b_, c_ = axpc.hist(counts[pos,1].astype(float), color = 'navy', density = True, bins = bins, cumulative = 1, alpha = 0.5)
print(a_)
print(b_)
if '--split_sets' in sys.argv:
    axpc.set_title('Positive correlated')
axpc.spines['top'].set_visible(False)
axpc.spines['right'].set_visible(False)
axpc.set_xlabel('Number drivers')
axpc.set_ylabel('Cumulative Number genes')
axpc.set_xlim([0,20])
axpc.grid()

if '--split_sets' in sys.argv:
    axnc = figc.add_subplot(122)
    axnc.hist(counts[neg,1].astype(float), color = 'navy', density = True, bins = bins, cumulative = 1, alpha = 0.5)
    axnc.set_title('Negative correlated')
    axnc.spines['top'].set_visible(False)
    axnc.spines['right'].set_visible(False)
    axnc.set_xlabel('Number drivers')
    axnc.set_ylabel('Cumulative Number genes')
    axnc.set_xlim([0,20])
    axnc.grid()


fig.savefig(os.path.splitext(sys.argv[1])[0]+'_distribution.jpg', dpi = 250, bbox_inches = 'tight')
figc.savefig(os.path.splitext(sys.argv[1])[0]+'_cumdistribution.jpg', dpi = 250, bbox_inches = 'tight')
print(os.path.splitext(sys.argv[1])[0]+'_distribution.jpg')

