# Plot gradient plot for GSTM3 with drivers and other variants
#python3 ../plot_attribution_alongsequence.py ENSG00000134202 0 ../ism_res/ENSG00000134202_ism_attributions.txt ../variant_info_100k/ENSG00000134202_frequency.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt ../ism_res/ENSG00000134202_ism_attributions_driversfw_types.txt --savefig

import numpy as np
import matplotlib.pyplot as plt
import logomaker as lm
import sys, os
import pandas as pd

# Gene name
gene = sys.argv[1]
# Plot array between position [window:-window]
window = int(sys.argv[2]) 
# Attributions and position of SNVs
varatts = np.genfromtxt(sys.argv[3], dtype = float)
# Population frequency of SNVs
popfreq = np.genfromtxt(sys.argv[4], dtype = float)
# align driver attributions and population frequency of snp
varatts, popfreq = varatts[np.argsort(varatts[:,0])[np.isin(np.sort(varatts[:,0]),popfreq[:,0])]], popfreq[np.argsort(popfreq[:,0])[np.isin(np.sort(popfreq[:,0]),varatts[:,0])]]

# normalize by maxmum attribution of any driver
maxvaratt = np.amax(np.absolute(varatts[:,-1]))
varatts[:,-1]/maxvaratt

# List of drivers for the gene
driverlist = open(sys.argv[5], 'r').readlines()
driverlist = np.array([line.strip().split() for line in driverlist])

# Drivertype for the gene
drivertype = open(sys.argv[6], 'r').readlines() # to select and assign type
drivertype = np.array([line.strip().split() for line in drivertype])
# align driverlist with driver type
driverlist = driverlist[np.argsort(driverlist[:,0])[np.isin(np.sort(driverlist[:,0]), drivertype[:,0])]]
drivertype = drivertype[np.argsort(drivertype[:,0])[np.isin(np.sort(drivertype[:,0]), driverlist[:,0])]]

loci = driverlist[:,0].astype(int)
attlen = driverlist[:,1].astype(float)
attlen = attlen/maxvaratt

# remove drivers from varatts
varatts = varatts[~np.isin(varatts[:,0],loci)]
drfreq = popfreq[np.isin(popfreq[:,0],loci),1]
popfreq = popfreq[~np.isin(popfreq[:,0],loci)]


allctypes = drivertype[:,-2:].astype(float)
bitypes = np.zeros((len(drivertype),2))
for c, ctypes in enumerate(allctypes):
    if (ctypes[0]<0 and ctypes[1]<0) or (ctypes[0]>0 and ctypes[1]>0):
        bitypes[c,0] = 1
    elif ctypes[0]>0 and ctypes[1]<0:
        bitypes[c,1] = 1
    elif ctypes[0]<0 and ctypes[1]>0:
        bitypes[c,1] = 1

tsslist = np.genfromtxt('../geneTSS.txt', dtype = str)

nts = list('ACGT')
ntsar = np.array(nts)


# Load gradient and one-hot encoded sequence for gene
grad = np.load('../gradient_tensors/'+gene+'_complete_grad_at_ref.npy')[0]
genehot = np.load('../ref_seqs/' +gene+'.npy').T
genehot = genehot[:, [0,3,2,1]]

lengene = np.shape(genehot)[0]

# determine the center of the gradient and the sequence file
centergrad = int(np.shape(grad)[0]/2) -1
centergene = int(np.shape(genehot)[0]/2)-1
# adjust the size of the gradient to the size of the sequence
offset = centergrad - centergene
grad = grad[offset:offset+lengene]

# find tss for the gene
tss = int(tsslist[list(tsslist[:,0]).index(gene),1])
# determine locations of drivers in the array
varintens = centergene + loci - tss
# determine locations of all other variants in the array
locvars = centergene + varatts[:,0].astype(int) - tss

grad /= np.amax(np.absolute(grad))

grad, genehot, varintens, centergene, locvars = grad[window:lengene-window], genehot[window:lengene-window], varintens - window, centergene - window, locvars - window
# remove variants that are located outside the sequence, ie. smaller than 0
varatts, popfreq, locvars = varatts[locvars>0], popfreq[locvars>0], locvars[locvars>0]

# Compute mean attribution for changing ref base to any other base
attribution = grad[genehot==1] - (np.sum(grad,axis = 1)-grad[genehot==1])/3


# compute mean and std with sliding windows
stepsize = 32
wdsize = 128
stdattribution = np.array([np.std(attribution[i:i+wdsize]) for i in range(0,len(attribution)-wdsize,stepsize)])
meanattribution = np.array([np.mean(attribution[i:i+wdsize]) for i in range(0,len(attribution)-wdsize,stepsize)])
stdx = np.arange(stepsize/2,len(attribution)-wdsize+stepsize/2,stepsize)

fig = plt.figure(figsize = (15,1.5))
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.fill_between(np.arange(len(attribution)), attribution, color = 'grey', label = 'Gradient attribution')
ax.plot([centergene, centergene],[np.amin(attribution),np.amax(attribution)], color = 'goldenrod', ls = ':')
if '--savefig' not in sys.argv:
    #ax.plot(np.arange(len(attribution)), attribution, color = 'grey', marker = '.')
    ax.scatter(np.arange(len(attribution)), attribution, color = 'grey', marker = '.')
    ax.plot([0,len(attribution)], [0.,0], color = 'k')
    ax.plot(stdx, stdattribution, color = 'navy', lw = 0.5)
    ax.plot(stdx, meanattribution, color = 'orange', lw = 0.5)
    ax.scatter(locvars, varatts[:,1], c = popfreq[:,1], vmin = 0, vmax = 1, cmap = 'Blues', edgecolor='silver')
else:
    ax.scatter(locvars, varatts[:,1], alpha = popfreq[:,1], vmin = 0, vmax = 1, color = 'b', edgecolor='none', marker = '.', label = 'SNPs')

maxdriver = np.argmax(driverlist[:,-3].astype(float))
for v, lo in enumerate(varintens):
    bounds = [0, attlen[v]]
    if types[v,0] == 1:
        marker = 's'
    elif types[v,1] == 1:
        marker = '^'
    elif types[v,2] == 1:
        marker = 'v'
    colr = 'red'
    label = None
    if v == maxdriver:
        colr = 'magenta'
        label = 'Main driver'
#    ax.plot([lo, lo],bounds, color = colr)
    ax.scatter([lo],[bounds[1]], cmap = 'Blues', vmin = 0, vmax = 1, c = drfreq[v], edgecolor = colr, marker = marker, label = label)

ax.set_xticks([centergene-40000,centergene,centergene+40000])
ax.set_xticklabels([-40000,'TSS',40000], rotation = 0)
ax.grid(color = 'grey', axis = 'y')
ax.set_ylabel('Attribution')
ax.set_xlim([0,len(attribution)])
#ax.legend()

if '--savefig' in sys.argv:
    fig.savefig(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseq.jpg', dpi = 200, bbox_inches = 'tight')
    print(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseq.jpg')
   
    ax.set_xticks([centergene-20000,centergene,centergene+20000])
    ax.set_xticklabels([-20000,'TSS',20000],rotation = 0)
    ax.set_xlim([centergene-30000, centergene+30000])
    fig.savefig(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseqcenter.jpg', dpi = 200, bbox_inches = 'tight')


    #ax.plot([0,len(attribution)], [0.,0], color = 'k')
    ax.set_xticks(np.append([centergene],varintens))
    ax.set_xticklabels(np.append(['TSS'],loci.astype(str)), rotation = 90)
    ax.scatter(np.arange(len(attribution)), attribution, color = 'grey', marker = '.')
    #ax.plot([0,len(attribution)], [0.,0], color = 'k')
    #ax.plot(stdx, stdattribution, color = 'blue', lw = 1)
    #ax.plot(stdx, meanattribution, color = 'orange', lw = 1)
    ax.scatter(locvars, varatts[:,1], c = popfreq[:,1], vmin = 0, vmax = 1, cmap = 'Blues', edgecolor='k')
    for v, lo in enumerate(varintens):
        bounds = [0, attlen[v]]
        if types[v,0] == 1:
            marker = 's'
        elif types[v,1] == 1:
            marker = '^'
        elif types[v,2] == 1:
            marker = 'v'
        colr = 'red'
        if v == maxdriver:
            colr = 'magenta'
        #ax.plot([lo, lo],bounds, color = colr)
        ax.scatter([lo],[bounds[1]], cmap = 'Blues', vmin = 0, vmax = 1, c = drfreq[v], edgecolor = colr, marker = marker)

    ax.set_xlim([np.amin(varintens)-100, np.amax(varintens)+100])
    fig.savefig(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseqclose.jpg', dpi = 200, bbox_inches = 'tight')
    print(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseqclose.jpg')
    
    ax.set_xticks([centergene-100,centergene,centergene+100])
    ax.set_xticklabels([-100,'TSS',100],rotation = 0)
    ax.set_xlim([centergene-150, centergene+150])
    fig.savefig(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseqtss.jpg', dpi = 200, bbox_inches = 'tight')
    print(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'_ingeneseqtss.jpg')
else:
    plt.show()



