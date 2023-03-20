# Read in gradient attributions along the entire sequence and cluster them on normalized standard deviation
# python3 ../cluster_grad_attributions.py ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list ../Prediction_correlationsCageAdultBrain.tx

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from sklearn.cluster import AgglomerativeClustering

# List of genes
genes = np.genfromtxt(sys.argv[1], dtype = str)
# Correlation of genes with observed across individuals
genecorr = np.genfromtxt(sys.argv[2], dtype = str)

# sort and select geneset
genes = np.sort(genes)
genecorr = genecorr[np.argsort(genecorr[:,0])[np.isin(np.sort(genecorr[:,0]), genes)]]
genecorr = genecorr[:,1].astype(float) > 0

nts = list('ACGT')
ntsar = np.array(nts)
window = 128 # window size to compute std

allstd = []
allatt = []
# read in gradient files and one-hot sequence encoding files and compute std and mena
for g, gene in enumerate(genes):
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

    attribution = np.sum(grad,axis = 1)

    attribution[np.sum(genehot,axis = 1) > 0] = grad[genehot==1] - (attribution[np.sum(genehot,axis = 1) > 0]-grad[genehot==1])/3
    attribution = attribution/np.std(attribution)
    allatt.append(attribution)
    stds = []
    for i in range(0, lengene-int(window/2), window):
        stds.append(np.std(attribution[i:i+window]))
    allstd.append(stds)
    
allstd = np.array(allstd)
allatt = np.array(allatt)

# Perform clustering and visualization of clusters
minsize = 0 # minsize could be adjusted to only show clusters with at least minsize genes
for nc in [2,5,10,20,40]: # Cluster data into 2,5,10,20,40 clusters with complete linkage
    clustering = AgglomerativeClustering(n_clusters = nc, affinity = 'euclidean', linkage = 'complete', distance_threshold= None).fit(allstd)
    clusters = clustering.labels_
    unclust, unclustn = np.unique(clusters, return_counts = True)
    unclust, unclustn = unclust[np.argsort(unclustn)], np.sort(unclustn)
    print(len(unclust), int(np.sum(unclustn>minsize)))
    height = 0.8/int(np.sum(unclustn>minsize))
    fig = plt.figure(figsize=(15,1*len(unclust[unclustn>minsize])), dpi = 100)
    axs = []
    axbs = []
    meanstds=[]
    maxattrib = []
    maxnum = []
    for i,c in enumerate(unclust[unclustn>minsize]):
        ax = fig.add_subplot(int(np.sum(unclustn>minsize)),1,i+1)
        ax.set_position([0.1, 0.9 - (i+0.95)*height, 0.8, height * 0.9])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        cmean = np.mean(np.absolute(allatt[clusters == c]), axis = 0)
        smean = np.mean(allstd[clusters == c], axis = 0)
        maxattrib.append(np.amax(cmean))
        meanstds.append(smean)
        ax.fill_between(np.arange(lengene), cmean, color = 'grey', alpha = 0.5)
        ax.plot(np.arange(window/2, lengene, window), smean)
        ax.plot([centergene, centergene],[np.amin(cmean), np.amax(cmean)], color= 'goldenrod', ls = '--')
        ax.plot([centergene-192, centergene-192],[np.amin(cmean), np.amax(cmean)/2], color= 'darkgoldenrod', ls = '--')
        ax.plot([centergene+192, centergene+192],[np.amin(cmean), np.amax(cmean)/2], color= 'darkgoldenrod', ls = '--')
        if i == int(np.sum(unclustn>minsize))-1:
            ax.set_xticks([centergene-20000, centergene, centergene +20000])
            ax.set_xticklabels(['-20,000', 'TSS', '20,000'], rotation = 60)
        else:
            ax.set_xticks([centergene-20000, centergene, centergene +20000])
            ax.tick_params(bottom = True, labelbottom = False)
        ax.set_xlim([0,lengene])
        axs.append(ax)
        axb = fig.add_subplot(int(np.sum(unclustn>minsize)),10,i+1)
        axb.set_position([0.91, 0.9 - (i+0.95)*height, 0.03, height * 0.9])
        axb.spines['top'].set_visible(False)
        axb.spines['left'].set_visible(False)
        axb.bar([0,1],[np.sum(genecorr[clusters == c]), np.sum(~genecorr[clusters == c])], width = 0.6, color = 'navy')
        axb.set_xticks([0,1])
        axb.set_xticklabels(['+', '-'])
        axb.tick_params(left = False, labelleft = False, right = True, labelright = True)
        maxnum.append(max(np.sum(~genecorr[clusters == c]), np.sum(genecorr[clusters ==c])))
        axb.set_yticks([maxnum[-1]])
        axbs.append(axb)
    print('ylim', np.median(maxattrib))
    for i, ax in enumerate(axs):
        ax.set_ylim([0.05, np.median(maxattrib)])
        axbs[i].set_ylim([0,np.amax(maxnum)])
    fig.savefig('Gradientstd128bp_clusters'+str(nc)+'.jpg', dpi = 200, bbox_inches = 'tight')
    print('Gradientstd128bp_clusters'+str(nc)+'.jpg')
    meanstds = np.array(meanstds)
    for m, mcut in enumerate([2,4,6]):
        mask = np.where(np.sum(meanstds>mcut,axis = 0)>0)[0]
        print(mask[0]*window-centergene, mask[-1]*window-centergene)
        for i, ax in enumerate(axs):
            ax.set_xlim([mask[0]*window, mask[-1]*window])
        fig.savefig('Gradientstd128bp_clusters'+str(nc)+'_zoomcut'+str(mcut)+'.jpg', dpi = 200, bbox_inches = 'tight')
        print('Gradientstd128bp_clusters'+str(nc)+'_zoomcut'+str(mcut)+'.jpg')
    for r in [10000, 5000, 2000, 1000]:
        for i, ax in enumerate(axs):
            ax.set_xlim([centergene-r, centergene+r])
            if i == int(np.sum(unclustn>minsize))-1:
                ax.set_xticks([min(centergene-r+500, centergene-r+1000), centergene, max(centergene +r-1000,centergene +r-500)])
                
                ax.set_xticklabels([str(int(min(centergene-r+500, centergene-r+1000)-centergene)), 'TSS', str(int(max(centergene +r-1000,centergene +r-500)-centergene))], rotation = 60)
            else:
                ax.set_xticks([centergene-r+1000, centergene, centergene+r-1000])
        fig.savefig('Gradientstd128bp_clusters'+str(nc)+'_zoomin'+str(r)+'.jpg', dpi = 200, bbox_inches = 'tight')
        print('Gradientstd128bp_clusters'+str(nc)+'_zoomin'+str(r)+'.jpg')



