# Generate scatter plot of computed pvalues versus the predicted correlation between observed and predicted expression values across individuals
# run 
#python3 scatter_pvalue_vs_correlationprediction.py --correlations Prediction_correlationsCageAdultBrain_Allstats.txt --alternativep Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --markersize 4 --printset 0

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import gaussian_kde, pearsonr

# Read correlations, y-axis values
pred = np.genfromtxt(sys.argv[1], dtype = str)

# Read p-values, x-axis values
pvalue = np.genfromtxt(sys.argv[2], dtype = str)
outname = os.path.splitext(sys.argv[2])[0]+'_to_PearsonR_scatter'

# Sort predicted correlaitions to match pvalues
prnames, prsort = np.unique(pred[:,0], return_index = True)
prsort = prsort[np.isin(prnames, pvalue[:,0])]
pred = pred[prsort]
# Sort pvalues to match predictions
pvsort = np.argsort(pvalue[:,0])[np.isin(np.sort(pvalue[:,0]), pred[:,0])]
pvalue = pvalue[pvsort]
# check if sorting was successful
print(np.array_equal(pred[:,0], pvalue[:,0]))

# select column with the data, use the mean predicted expression as the coloring
predicted, colors = pred[:,1].astype(float), np.abs(pred[:,-2].astype(float)) #-pred[:,-2].astype(float)) 
# select p-value as x-axis data
measured = pvalue[:,-2].astype(float)
# determine cut-off from corrected pvalue with FDR <0.05
cut = np.amin(np.abs(measured[pvalue[:,-1].astype(float)>=-np.log10(0.05)]))
print(cut, 10**-cut)


mask = ~np.isnan(predicted) * ~np.isnan(colors) * ~np.isnan(measured)
measured, predicted, colors, pred = measured[mask], predicted[mask], colors[mask], pred[mask]

print(len(measured))

for s in np.argsort(predicted)[::-1][:20]:
    print(pred[s,0], measured[s], predicted[s])
print('\n')
for s in np.argsort(predicted)[:20]:
    print(pred[s,0], measured[s], predicted[s])

cmap = 'inferno'
vmin, vmax = 0,  int(np.amax(np.abs(colors)))
print(np.amin(colors))
if np.amin(colors) < 0:
    vmin = -vmax

# Choose to color a set of genes darker than the unselected set
if '--genesetcolor' in sys.argv:
    colors = np.zeros(len(measured))+0.2
    geneset = np.genfromtxt(sys.argv[sys.argv.index('--genesetcolor')+1], dtype = str)
    colors[np.isin(pred[:,0], geneset)] = 1.15
    outname += os.path.splitext(os.path.split(sys.argv[sys.argv.index('--genesetcolor')+1])[1])[0]
    vmin, vmax = 0,  int(np.ceil(np.amax(colors)))
# Or color scatters based on a continuous value
elif '--colors' in sys.argv:
    colors = np.genfromtxt(sys.argv[sys.argv.index('--colors')+1], dtype = str)
    colors = colors[np.argsort(colors[:,0])[np.isin(np.sort(colors[:,0]), pred[:,0])]]
    print(np.array_equal(colors[:,0], pred[:,0]))
    print(len(colors), len(pred))
    colors = np.nan_to_num(colors[:, int(sys.argv[sys.argv.index('--colors')+2])].astype(float))
    outname += os.path.splitext(os.path.split(sys.argv[sys.argv.index('--colors')+1])[1])[0]+sys.argv[sys.argv.index('--colors')+2]
    if np.amin(colors) < 0:
        vmin, vmax = -np.amax(np.abs(colors)), np.amax(np.abs(colors))
        cmap = 'BrBG'
    else:
        vmin, vmax = 0, np.amax(np.abs(colors))
    # Set amin and amax based on 99 quartile
    if '--quartilecolor' in sys.argv:
        if np.amin(colors) < 0:
            vmin, vmax = -np.percentile(np.abs(colors),99), np.percentile(np.abs(colors), 99)
            cmap = 'BrBG'
        else:
            vmin, vmax = 0, np.percentile(np.abs(colors),99)
    print(vmin,vmax)




# Adjust different pyplot features
if '--colormap' in sys.argv:
    cmap = sys.argv[sys.argv.index('--colormap')+1]

if '--markersize' in sys.argv:
    plt.rcParams['lines.markersize'] = float(sys.argv[sys.argv.index('--markersize')+1])
lw = 0
if '--linewidth' in sys.argv:
    lw = 0.3
alpha = 1
if '--alpha' in sys.argv:
    alpha = float(sys.argv[sys.argv.index('--alpha')+1])

fig = plt.figure(figsize = (4,4.), dpi = 200)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
sort = np.argsort(np.abs(colors))
print('Mean', np.mean(predicted[sort]))
print('PosSig', int(np.sum((predicted>0) * (measured>cut))), 'NegSig', int(np.sum((predicted<0) * (measured>cut))))
print('PosSig0.1', int(np.sum((predicted>0.1) * (measured>cut))), 'NegSig', int(np.sum((predicted<-0.1) * (measured>cut))))
print('PosSig0.2', int(np.sum((predicted>0.2) * (measured>cut))), 'NegSig', int(np.sum((predicted<-0.2) * (measured>cut))))

a = ax.scatter(measured[sort], predicted[sort], cmap = cmap, vmin = vmin, vmax = vmax, c=colors[sort], alpha = alpha, marker = 'o', lw = lw, edgecolor ='silver' )
ax.plot([cut,cut],[np.amin(predicted), np.amax(predicted)], color = 'r', ls = '--')
ax.plot([np.amin(measured), np.amax(measured)], [0, 0], color = 'grey', ls = '--')
ax.set_xlabel('Log10(p-value)')
ax.set_ylabel('Pearson R')
if '--genesetcolor' not in sys.argv:
    fig.colorbar(a, pad = -0.1, fraction = 0.09, shrink = 0.25, aspect = 2, anchor = (0.,0.99), ticks = [0, int(vmax)], label = None)
print(outname)
fig.savefig(outname + '.jpg', dpi = 450, bbox_inches = 'tight')

# Print out set of genes that are above cut and also have abs(y-axis) > pcut
if '--printset' in sys.argv:
    pcut = float(sys.argv[sys.argv.index('--printset')+1])
    for s in np.where((measured>cut) * (np.abs(predicted) > pcut))[0]:
        print(pred[s,0], measured[s], predicted[s], colors[s])


