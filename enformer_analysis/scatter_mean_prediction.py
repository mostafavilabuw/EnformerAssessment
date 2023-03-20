# Generate a scatter plot between predicted and observed expression values

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import gaussian_kde, pearsonr

pred = np.genfromtxt('MeanGeXPredFineTuned.txt', dtype = str, skip_header = 1)
genes, measured, predicted = pred[:,0], pred[:,1].astype(float), pred[:,2].astype(float)

outname = 'MeanExpression_to_referencePrediction_scatter.jpg'
if '--testset' in sys.argv:
    tset = np.load(sys.argv[sys.argv.index('--testset')+1])
    mask = np.isin(genes, tset)
    genes, measured, predicted = genes[mask], measured[mask], predicted[mask]
    outname = 'MeanExpression_to_referencePrediction_scatter_testset.jpg'

colors = gaussian_kde(np.array([measured, predicted]))(np.array([measured, predicted]))

fig = plt.figure(figsize = (4.,3.), dpi = 200)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
sort = np.argsort(colors)
ax.scatter(measured[sort], predicted[sort], cmap = 'cividis', c=colors[sort], alpha = 0.5, label = 'R='+str(round(pearsonr(measured, predicted)[0],2)))
ax.set_xlim([2, np.amax(measured)*1.01])
ax.legend()
ax.set_xlabel('Observered mean expression')
ax.set_ylabel('Predicted reference expression')
fig.savefig(outname, dpi = 450, bbox_inches = 'tight')



