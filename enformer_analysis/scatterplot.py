# script to plot to values against each other and color color in different ways
# python3 scatterplot.py Prediction_correlationsCageAdultBrain_Allstats.txt ism_res/Refpred.txt "CAGE,brain,adult,MeanSumlog10+-2indv" "CAGE,brain,adult,log10sum+-1ref" --columns -2 -1 --density --alpha 0.5 --label --filternan --linewidth 0.1 --log10y

import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr 
from scipy.spatial.distance import cdist

def read(file, column = -1):
    indv, val = [], []
    for l, line in enumerate(open(file, 'r')):
        if line[0] != '#':
            line = line.strip().split()
            indv.append(line[0])
            val.append(float(line[column]))
    indv, val = np.array(indv), np.array(val)
    _, sort = np.unique(indv, return_index = True)
    indv, val = indv[sort], val[sort]
    return indv, val

# Columns that will be selected from the txt files
column0, column1 = -1, -1
if '--columns' in sys.argv:
    column0, column1 = int(sys.argv[sys.argv.index('--columns')+1]), int(sys.argv[sys.argv.index('--columns')+2])
# Read two text files
print(sys.argv[1], column0)
indv1, val1 = read(sys.argv[1], column0)
indv2, val2 = read(sys.argv[2], column1)

# Determinne labels of x and y-axis
xname = sys.argv[3]
yname = sys.argv[4]

outname = os.path.splitext(sys.argv[1])[0]+'_vs_'+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+'_scatter'
if column0 != -1 or column1 != -1:
    outname += str(column0)+'-'+str(column1)

# sort and align values to each other
sort1 = np.isin(indv1, indv2)
sort2 = np.isin(indv2, indv1)

indv1, val1 = indv1[sort1], val1[sort1]
indv2, val2 = indv2[sort2], val2[sort2]

# check if sorting was successful
print(np.array_equal(indv1, indv2), len(indv1), len(indv2))

# remove nans
if '--filternan' in sys.argv:
    mask = ~np.isnan(val1) * ~np.isnan(val2)
    indv1, indv2, val1, val2 = indv1[mask], indv2[mask], val1[mask], val2[mask]

if '--minx' in sys.argv:
    mask = val1 > float(sys.argv[sys.argv.index('--minx')+1])
    indv1, indv2, val1, val2 = indv1[mask], indv2[mask], val1[mask], val2[mask]

if '--miny' in sys.argv:
    mask = val2 > float(sys.argv[sys.argv.index('--miny')+1])
    indv1, indv2, val1, val2 = indv1[mask], indv2[mask], val1[mask], val2[mask]


if '--absy' in sys.argv:
    val2 = np.abs(val2)
if '--absx' in sys.argv:
    val1 = np.abs(val1)
if '--log10y' in sys.argv:
    val2 = np.log10(val2+1)
if '--log10x' in sys.argv:
    val1 = np.log10(val1+1)


if '--percentilex' in sys.argv:
    perc = float(sys.argv[sys.argv.index('--percentilex')+1])
    pmin, pmax = np.percentile(val1, [100.-perc, perc])
    mask = (val1 < pmax) * (val1 > pmin)
    indv1, indv2, val1, val2 = indv1[mask], indv2[mask], val1[mask], val2[mask]

if '--percentiley' in sys.argv:
    perc = float(sys.argv[sys.argv.index('--percentiley')+1])
    pmin, pmax = np.percentile(val2, [100.-perc, perc])
    mask = (val2 < pmax) * (val2 > pmin)
    indv1, indv2, val1, val2 = indv1[mask], indv2[mask], val1[mask], val2[mask]


colors = np.ones(len(indv1))*0.6
edgecolors = np.chararray(len(indv1), itemsize = 20, unicode = True)
edgecolors[:] = 'grey'
markers = np.array(['o' for i in range(len(indv1))])

# assign marker types based on drivers or other features
if '--assigntype' in sys.argv:
    drivertype = open(sys.argv[sys.argv.index('--assigntype')+1], 'r').readlines() # to select and assign type
    if len(drivertype) > 0:
        drivertype = np.array([line.strip().split() for line in drivertype])
        allctypes = drivertype[:,-2:].astype(float)
        drivertype = drivertype[:,0] #.astype(int)
        # types:
        # 0: correct (ar<0&eq<0 or ar>0&eq>0)
        # 1: false positive (ar>0&eq<0)
        # 2: false negative (ar<0&eq>0)
        for c, ctypes in enumerate(allctypes):
            if (ctypes[0]<0 and ctypes[1]<0) or (ctypes[0]>0 and ctypes[1]>0):
                edgecolors[list(indv1).index(drivertype[c])] = 'r'
                markers[list(indv1).index(drivertype[c])] = 's'
            elif ctypes[0]>0 and ctypes[1]<0:
                edgecolors[list(indv1).index(drivertype[c])] = 'r'
                markers[list(indv1).index(drivertype[c])] = '^'
            elif ctypes[0]<0 and ctypes[1]>0:
                edgecolors[list(indv1).index(drivertype[c])] = 'r'
                markers[list(indv1).index(drivertype[c])] = 'v'

# assign main driver
if '--assignmain' in sys.argv:
    maindriver = open(sys.argv[sys.argv.index('--assignmain')+1], 'r').readlines() # to select and assign type
    if len(drivertype) > 0:
        maindriver = np.array([line.strip().split() for line in maindriver])
        edgecolors[list(indv1).index(maindriver[np.argmax(maindriver[:,-3].astype(float)),0])] = 'magenta'

vmin, vmax = 0, 1
cmap = 'Blues'
# assign color to dots
if '--colors' in sys.argv:
    colors = np.genfromtxt(sys.argv[sys.argv.index('--colors')+1], dtype = str)
    colors = colors[np.argsort(colors[:,0])[np.isin(np.sort(colors[:,0]),indv1)]]
    colors = colors[:,int(sys.argv[sys.argv.index('--colors')+2])].astype(float)
    if np.amin(colors) < 0:
        vmin, vmax = -np.amax(np.absolute(colors)), np.amax(np.absolute(colors))
        cmap = 'RdBu_r'
# color based on density of dots
if '--density' in sys.argv:
    from scipy.stats import gaussian_kde
    colors = gaussian_kde(np.array([val1, val2]))(np.array([val1, val2]))
    vmin, vmax = np.amin(colors), np.amax(colors)
    cmap = 'viridis'
# If clusters are loaded, all dots in a cluster will be connected by a line
if '--connect_clusters' in sys.argv:
    clusters = np.genfromtxt(sys.argv[sys.argv.index('--connect_clusters')+1], dtype = str)
    clusters = clusters[np.argsort(clusters[:,0])[np.isin(np.sort(clusters[:,0]),indv1)]]
    
fig = plt.figure(figsize = (3.8,3.5), dpi = 200)
ax =  fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

if np.amin(val1) < 0:
    ax.plot([0,0],[np.amin(val2), np.amax(val2)], color = 'k', ls = '--')
if np.amin(val2) < 0:
    ax.plot([np.amin(val1), np.amax(val1)],[0,0], color = 'k', ls = '--')

if '--connect_clusters' in sys.argv:
    ucluster, nucluster = np.unique(clusters, return_counts = True)
    ucluster = ucluster[nucluster > 1]
    for c in ucluster:
        allgood = True
        mask = np.where(clusters == c)[0]
        if '--assigntype' in sys.argv:
            if (edgecolors[mask] != 'grey').any():
                allgood = True
            else:
                allgood = False
        if len(mask) > 2 and allgood:
            distmat = cdist(np.array([val1[mask], val2[mask]]).T,np.array([val1[mask], val2[mask]]).T ,'euclidean')
            maxdist = np.amax(distmat)+1
            np.fill_diagonal(distmat, maxdist)
            count = np.array(np.where(distmat == np.amin(distmat))).T[0]
            hascon = np.zeros(len(distmat))
            ax.plot(val1[mask][count], val2[mask][count], color = 'royalblue', lw = 0.5)
            distmat[count] = distmat[count[::-1]] = maxdist
            hascon[count] += 1
            while True:
                conn = np.array(np.where(distmat == np.amin(distmat[count][hascon[count]<2]))).T
                connmask = np.sum(np.isin(conn, count),axis =1) == 1
                conn = conn[connmask]
                conn = conn[np.sum(np.isin(conn,np.where(hascon < 2)[0]),axis = 1)==2]
                conn = conn[0]
                ax.plot(val1[mask][conn], val2[mask][conn], color = 'royalblue', lw = 0.5)
                distmat[conn] = distmat[conn[::-1]] = maxdist
                hascon[conn] += 1
                count = np.unique(np.append(count, conn))
                if len(count) == len(mask):
                    break
        elif allgood:
            ax.plot(val1[mask], val2[mask], color = 'royalblue', lw = 0.5)


lw = 0.5
if '--linewidth' in sys.argv:
    lw = float(sys.argv[sys.argv.index('--linewidth')+1])

alpha = 1.
if '--alpha' in sys.argv:
    alpha = float(sys.argv[sys.argv.index('--alpha')+1])

label = None
if '--label' in sys.argv:
    pears, pval = pearsonr(val1, val2)
    label = 'R='+str(round(pears,3))+'\np='+str(round(pval,4))
marksort, marksortn = np.unique(markers, return_counts = True)
marksort = marksort[np.argsort(-marksortn)]
for mark in marksort:
    mask = np.where(markers == mark)[0]
    mask = mask[np.argsort(np.absolute(colors[mask]))]
    ax.scatter(val1[mask], val2[mask], cmap = cmap, vmin = vmin, vmax = vmax, c = colors[mask], alpha = alpha, lw = lw, edgecolor = list(edgecolors[mask]), marker = mark, label = label)

ax.set_xlabel(xname)
ax.set_ylabel(yname)
if label is not None:
    ax.legend()

# Print stats from every quadrant
if '--print_quadrant1' in sys.argv:
    cuts = sys.argv[sys.argv.index('--print_quadrant1')+1].split(',')
    mask = (val1<float(cuts[0])) * (val2>float(cuts[1]))
    for i in np.where(mask)[0]:
        print(indv1[i], val1[i], val2[i])
    
if '--print_quadrant2' in sys.argv:
    cuts = sys.argv[sys.argv.index('--print_quadrant2')+1].split(',')
    mask = (val1>float(cuts[0])) * (val2>float(cuts[1]))
    for i in np.where(mask)[0]:
        print(indv1[i], val1[i], val2[i])

if '--print_quadrant3' in sys.argv:
    cuts = sys.argv[sys.argv.index('--print_quadrant3')+1].split(',')
    mask = (val1>float(cuts[0])) * (val2<float(cuts[1]))
    for i in np.where(mask)[0]:
        print(indv1[i], val1[i], val2[i])

if '--print_quadrant4' in sys.argv:
    cuts = sys.argv[sys.argv.index('--print_quadrant4')+1].split(',')
    mask = (val1<float(cuts[0])) * (val2<float(cuts[1]))
    for i in np.where(mask)[0]:
        print(indv1[i], val1[i], val2[i])

if '--logx' in sys.argv:
    ax.set_xscale('log')
elif '--symlogx' in sys.argv:
    ax.set_xscale('symlog')
if '--logy' in sys.argv:
    ax.set_yscale('log')
elif '--symlogy' in sys.argv:
    ax.set_yscale('symlog')

if '--show' in sys.argv:
    plt.show()
else:
    fig.savefig(outname+'.jpg', dpi = 250, bbox_inches = 'tight')
    print(outname+'.jpg')






