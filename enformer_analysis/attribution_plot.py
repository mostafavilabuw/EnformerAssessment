# Plots attributions of SNVs along the genomic location
# Run:
# python3 ../attribution_plot.py ENSG00000013573_ism_attributions.txt ISM ../variant_info_100k/ENSG00000013573snp_info.txt --tss ../geneTSS.txt ENSG00000013573

import numpy as np
import sys, os
import matplotlib.pyplot as plt

# Read attributions
attributions = np.genfromtxt(sys.argv[1])
# Define ylabel
ylabel = sys.argv[2]

# Color attributions by other features such as population frequency
if '--colors' in sys.argv:
    colors = np.genfromtxt(sys.argv[sys.argv.index('--colors')+1])
    # sort colors to attributions
    sort = np.argsort(colors[:,0])[np.isin(np.sort(colors[:,0]), attributions[:,0])]
    colors = colors[sort, 1]

# Change vmin, vmax for colors if needed
vmin, vmax = 0, 1
if '--vmin' in sys.argv:
    vmin = int(sys.argv[sys.argv.index('--vmin')+1])
if '--vmax' in sys.argv:
    vmax = int(sys.argv[sys.argv.index('--vmax')+1])


figsize = (15,2) # None
fig = plt.figure(figsize = figsize)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
sort = np.argsort(colors)
a = ax.scatter(attributions[sort,0], attributions[sort,1], cmap = 'Blues', vmin = vmin, vmax = vmax, c = colors[sort], edgecolor = 'grey')
fig.colorbar(a, pad = 0.01, fraction = 0.09, shrink = 0.25, aspect = 2, anchor = (0.,0.9))

# Show markers with increased size and red edgecolor
if '--drivers' in sys.argv:
    dobj = open(sys.argv[sys.argv.index('--drivers')+1], 'r').readlines()
    dfile = [line.strip().split() for line in dobj]
    dfile = np.array(dfile,dtype =float)
    dloc = np.where(np.isin(attributions[:,0], dfile[:,0]))
    dloc = dloc[0]
    size0 = plt.rcParams['lines.markersize'] ** 2
    ax.scatter(attributions[dloc,0], attributions[dloc,1], s = size0 * (1+2*dfile[:,-3]), linewidths = 1.2, cmap = 'Blues', vmin = 0, vmax =1, c = colors[dloc], edgecolor = 'red')
    for t in dloc:
        print(attributions[t,0])
        ax.text(attributions[t,0], attributions[t,1], str(int(attributions[t,0])), ha = 'left', va = 'bottom')

else:
    dloc = np.argsort(-np.absolute(attributions[:,1]))[:6]
    # Add the name of the top attributions to the figure
    if '--name_top_attributions' in sys.argv:
        for t in dloc:
            print(attributions[t,0])
            ax.text(attributions[t,0], attributions[t,1], str(int(attributions[t,0])), ha = 'left')

# Mark the location of the TSS in the figure
if '--tss' in sys.argv:
    tss = np.genfromtxt(sys.argv[sys.argv.index('--tss')+1], dtype = str)
    tss = int(tss[list(tss[:,0]).index(sys.argv[sys.argv.index('--tss')+2]),1])
    ax.plot([tss, tss],[np.amin(attributions[:,1]), np.amax(attributions[:,1])], ls = '--', color = 'goldenrod', alpha = 0.8)
    ax.set_title(sys.argv[sys.argv.index('--tss')+2])
    xmin, xmax = np.amin(attributions[:,0]), np.amax(attributions[:,0])
    ax.set_xticks([tss-int((xmax-xmin)/40000)*10000, tss, tss+int((xmax-xmin)/40000)*10000])
    ax.set_xticklabels([-int((xmax-xmin)/40000)*10000, 'TSS', int((xmax-xmin)/40000)*10000])

ax.set_xlabel('Genomic location')
ax.set_ylabel(ylabel)
mina, maxa = np.amin(attributions[:,1]), np.amax(attributions[:,1])
dista = maxa - mina
ax.set_ylim([mina-0.1*dista, maxa + 0.1*dista])
ax.set_xlim([np.amin(attributions[:,0]), np.amax(attributions[:,0])])
fig.savefig(os.path.splitext(sys.argv[1])[0]+'.jpg', dpi = 250, bbox_inches = 'tight')
if '--tss' in sys.argv:
    zeroline = ax.plot([xmin, xmax], [0,0], color = 'grey', alpha = 0.5, zorder = -1)
    ax.set_xlim([tss-4000, tss+4000])
    ax.set_xticks([tss-2000, tss+2000])
    ax.set_xticklabels([-2000,2000])
    fig.savefig(os.path.splitext(sys.argv[1])[0]+'tsszoom.jpg', dpi = 350, bbox_inches = 'tight')

ax.set_xlim([np.amin(attributions[dloc,0])-2000, np.amax(attributions[dloc,0])+2000])
fig.savefig(os.path.splitext(sys.argv[1])[0]+'zoom.jpg', dpi = 350, bbox_inches = 'tight')

