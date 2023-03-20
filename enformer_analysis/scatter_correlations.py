# Plots a scatter plot between the correlation of two predictors to observed expression
#python3 scatter_correlations.py --colors Full_analysis/Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats.txt --colorscut 1.301 --vlim -1.5,1

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import gaussian_kde, pearsonr

# Read in correlation of Enformer and PrediXcan
pred = np.genfromtxt('Full_analysis/Prediction_correlationsCageAdultBrain_Allstats.txt', dtype = str)
predX = np.genfromtxt('PrediXcanCorrelationWithExpressionENSG.tsv', dtype = str)
predX = predX[predX[:,1]!= 'NaN']

# sort genes of both files so that rows match
a_, sortp = np.unique(pred[:,0], return_index = True)
pred = pred[sortp[np.isin(a_, predX[:,0])]]
a_, sortx = np.unique(predX[:,0], return_index = True)
predX = predX[sortx[np.isin(a_, pred[:,0])]]
# check sorting
print(np.array_equal(pred[:,0],predX[:,0]))

# Select rows that are plotted against each other
measured, predicted, colors = pred[:,3].astype(float), pred[:,1].astype(float), pred[:,4].astype(float)
mask = ~np.isnan(predicted)
measured, predicted, colors, pred, predX = measured[mask], predicted[mask], colors[mask], pred[mask], predX[mask]

# set cutoff as R>0.2 ~ pvalue 0.05
#cutoff = np.amin(np.absolute(predicted[pred[:,2].astype(float)<0.01]))
cutoff = 0.2
predX = predX[:,1].astype(float)

outname = 'PearsonR_PredixCan_Enformer_scatter.jpg'
outname2 = 'PearsonR_PredixCan_absEnformer_scatter.jpg'


colors = np.chararray(len(colors), itemsize = 20, unicode = True)
colors[:] = 'indigo'
# Choose colors 
if '--colors' in sys.argv:
    colorfile = sys.argv[sys.argv.index('--colors')+1]
    colors = np.genfromtxt(colorfile, dtype = str)
    ic_, csort = np.unique(colors[:,0], return_index = True)
    colors = colors[csort[np.isin(colors[csort,0], pred[:,0] )]]
    colors = colors[colors[:,-1]!='nan']
    
    outname = os.path.splitext(outname)[0] + os.path.splitext(os.path.split(colorfile)[1])[0] + '.jpg'
    outname2 = os.path.splitext(outname2)[0] + os.path.splitext(os.path.split(colorfile)[1])[0] + '.jpg'
    
    pnames, mask = np.unique(pred[:,0], return_index = True)
    mask = mask[np.isin(np.sort(pnames),colors[:,0])]
    measured, predicted, pred, predX = measured[mask], predicted[mask], pred[mask], predX[mask]
    
    print(np.array_equal(colors[:,0], pred[:,0]), len(colors), len(pred))
    colors = colors[:,-1].astype(float)
    # make color binary
    if '--colorscut' in sys.argv:
        colcut = float(sys.argv[sys.argv.index('--colorscut')+1])
        colors[np.absolute(colors)< colcut] = 0
        colors[colors>= colcut] = 1
        colors[colors<= -colcut] = -1
        print(np.sum(colors == 1), np.sum(colors == -1))
        print(np.sum((np.absolute(colors)==1)*(predicted >0) ))
    # cutoff high color values
    elif '--colorsmax' in sys.argv:
        colcut = float(sys.argv[sys.argv.index('--colorsmax')+1])
        colors[np.abs(colors)>= colcut] = colcut
    elif '--logcolors' in sys.argv:
        colors = np.sign(colors)*np.log(np.abs(colors)+1)

# generate colors from density 
elif '--density' in sys.argv:
    colors = gaussian_kde(np.array([measured, np.abs(predicted)]))(np.array([measured, predicted]))

# Set other matplotlib parameters for the scatter plot
vmin, vmax = -1, 1
if '--vlim' in sys.argv:
    vmin, vmax = sys.argv[sys.argv.index('--vlim')+1].split(',')
    vmin, vmax = float(vmin), float(vmax)

lw = 0
if '--linewidth' in sys.argv:
    lw = float(sys.argv[sys.argv.index('--linewidth')+1])

alpha = 0.7
if '--alpha' in sys.argv:
    alpha = float(sys.argv[sys.argv.index('--alpha')+1])


cmap = 'Purples'
if '--cmap' in sys.argv:
    cmap = sys.argv[sys.argv.index('--cmap')+1]

# print some stats
print('\nHighest Enformer')
for s in np.argsort(predicted)[::-1][:10]:
    print(pred[s,0], predX[s], predicted[s])

print('\nLowest Enformer')
for s in np.argsort(predicted)[:10]:
    print(pred[s,0], predX[s], predicted[s])

print('\nHighest PrediXcan')
for s in np.argsort(predX)[::-1][:20]:
    print(pred[s,0], predX[s], predicted[s])

print('\nMean', np.mean(predicted), np.mean(predX))
print('Above', cutoff, np.sum(predicted>cutoff), np.sum(predX>cutoff), 'out of', len(predicted))

bothmask = (np.abs(predicted)>cutoff) * (predX>cutoff)
onlyx = (np.abs(predicted)<cutoff) * (predX>cutoff)
nonex = (np.abs(predicted)<cutoff) * (predX<cutoff)

print('\nboth above', cutoff)
for s in np.where(bothmask)[0]:
    print(pred[s,0])#, predX[s], predicted[s]) 

print('\nOnly PredX above', cutoff)
for s in np.where(onlyx)[0]:
    print(pred[s,0])#, predX[s], predicted[s])

print('\nNone', cutoff)
for s in np.where(nonex)[0]:
    print(pred[s,0])#, predX[s], predicted[s])

fig = plt.figure(figsize = (4,4.), dpi = 200)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
sort = np.argsort(colors)
a = ax.scatter(predX[sort], predicted[sort], cmap = cmap, c = colors[sort], alpha = alpha, marker = 'o', lw = lw, s = 20, edgecolor = 'silver', vmin = vmin, vmax = vmax)
ax.scatter([np.mean(predX)],[np.mean(predicted)], marker = 'x', color = 'limegreen')
ax.plot([np.amin(predX), np.amax(predX)], [cutoff, cutoff], color = 'r', ls = '--')
ax.plot([np.amin(predX), np.amax(predX)], [-cutoff, -cutoff], color = 'r', ls = '--')
ax.plot([np.amin(predX), np.amax(predX)], [0, 0], color = 'grey', ls = '-')
ax.plot([cutoff, cutoff], [np.amin(predicted), np.amax(predicted)], color = 'red', ls = '--')
ax.plot([0,0], [np.amin(predicted), np.amax(predicted)], color = 'grey', ls = '-')
ax.set_xlabel('PrediXcan R')
ax.set_ylabel('Enformer R')
fig.savefig(outname, dpi = 450, bbox_inches = 'tight')

fig = plt.figure(figsize = (4,4.), dpi = 200)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
sort = np.argsort(colors)
a = ax.scatter(predX[sort], np.abs(predicted[sort]), cmap = cmap, c = colors[sort], alpha = alpha, marker = 'o', lw = lw, s = 20, edgecolor = 'silver', label = 'R='+str(round(pearsonr(predX[sort], np.abs(predicted[sort]))[0],2)), vmin = vmin, vmax = vmax)
ax.scatter([np.mean(predX)],[np.mean(np.abs(predicted))], marker = 'x', color = 'limegreen')
ax.plot([np.amin(predX), np.amax(predX)], [cutoff, cutoff], color = 'r', ls = '--')
ax.plot([np.amin(predX), np.amax(predX)], [0, 0], color = 'grey', ls = '-')
ax.plot([cutoff, cutoff], [0, np.amax(predicted)], color = 'red', ls = '--')
ax.plot([0, np.amax(predicted)], [0, np.amax(predicted)], color = 'red', ls = '-')
ax.set_xlabel('PrediXcan R')
ax.set_ylabel('Enformer abs(R)')
ax.legend()
#fig.colorbar(a, pad = 0., fraction = 0.09, shrink = 0.25, aspect = 2, anchor = (0.,0.9), ticks = [0, int(np.amax(colors))], label = None)
fig.savefig(outname2, dpi = 450, bbox_inches = 'tight')




