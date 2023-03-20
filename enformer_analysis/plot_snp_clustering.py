# Visualize the SNV distribution across individuals
# Can also plot a clustered/sorted heatmap for the distribution of SNVs, i.e. SNVs versus individuals
# python3 ../plot_snp_clustering.py ENSG00000013573snp_info.txt ENSG00000013573row_names.txt --select_individuals ../Enformer_predictions_individuals.txt --minfreq 0.1 --minsnps 10 --combine_snps 0.9 --combine_individuals 0.9 --savefig ENSG00000013573_individual_enformer_snp


import numpy as np
import sys, os
import matplotlib.pyplot as plt
from functools import reduce
import glob
from matplotlib import cm
from scipy.spatial.distance import pdist, cdist
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import cdist, pdist
from sklearn.cluster import AgglomerativeClustering as agc

# Function to plot a sorted heatmap
# X and Y axis can get additional attributes and can be sorted according to different algorithms
def plot_heatmap(heatmat, measurex = None, measurey = None, sortx = None, sorty = None, x_attributes = None, y_attributes = None, xattr_name = None, yattr_name = None, heatmapcolor = cm.BrBG_r, xatt_color = None, yatt_color = None, pwms = None, combine_cutx = 0., combine_cuty = 0., color_cutx = 0., color_cuty = 0., plot_value = False, vmin = None, vmax = None, grid = False, xdenline = None, ydenline = None, xlabel = None, ylabel = None, xticklabels = None, yticklabels  = None, dpi = 100, figname = None, maxsize = 20, ratio = (16,9)):
    
    # either provide similarity matrix as heatmap (measurex = None) or provide a similarity function from scipy.spatial.distance.pdist
    # If no measure is provided heatmap entries will be rescaled between 0,1 and a transformation function can retransform for xticklabels
    if measurex is not None:
        simatrixX = pdist(heatmat.T, metric = measurex)
    else:
        if np.shape(heatmat)[0] != np.shape(heatmat)[1]:
            print( 'heatmat not symmetric matrix: sortx set to None if given')
            sortx = None
        else:
            if np.any(np.abs(heatmat - heatmat.T) > 1e-8):
                print( 'heatmat not symmetric matrix: sortx set to None if given')
                sortx = None
        
        if sortx is not None:        
            # checks if similarity matrix or distance matrix
            issimilarity = np.all(np.amax(heatmat) == np.diag(heatmat))
            heatmax, heatmin = np.amax(heatmat), np.amin(heatmat)
            simatrixX = 1.-heatmat.T #int(issimilarity) - (2.*int(issimilarity)-1.) * (heatmat - heatmin)/(heatmax - heatmin)
            simatrixX = simatrixX[np.triu_indices(len(simatrixX),1)]
            
    if measurey is not None:
        simatrixY = pdist(heatmat, metric = measurey)
    else:
        if np.shape(heatmat)[0] != np.shape(heatmat)[1]:
            print( 'heatmat not symmetric matrix: sorty set to None if given')
            sorty = None
        else:
            if np.any(np.abs(heatmat - heatmat.T) > 1e-8):
                print( 'heatmat not symmetric matrix: sorty set to None if given')
                sorty = None
        if sorty is not None:        
            # checks if similarity matrix or distance matrix
            issimilarity = np.all(np.amax(heatmat) == np.diag(heatmat))
            heatmax, heatmin = np.amax(heatmat), np.amin(heatmat)
            simatrixY = 1.-heatmat #int(issimilarity) - (2.*int(issimilarity)-1.) * (heatmat - heatmin)/(heatmax - heatmin)
            simatrixY = simatrixY[np.triu_indices(len(simatrixY),1)]
            
    # Generate dendrogram for x and y
    #### NJ not yet included
    if sortx is not None:
        Zx = linkage(simatrixX, sortx)
        #if combine_cutx > 0:
            #Zx = reduce_z(Zx, combine_cutx)
    if sorty == 'maxdist':
        sortsize = np.argsort(heatmat[:,0] -heatmat[:,1] + np.amax(heatmat, axis = 1)*0.1)
        
    elif sorty is not None:
        Zy = linkage(simatrixY, sorty) 
        #if combine_cuty > 0:
            #Zy = reduce_z(Zy, combine_cuty)
    xextra = 0.
    if y_attributes is not None:
        xextra = np.shape(y_attributes)[1]*0.8
    yextra = 0.
    if x_attributes is not None:
        yextra = np.shape(x_attributes)[0]*0.8
    

    fig = plt.figure(figsize = (min(maxsize*ratio[0]/ratio[1],0.3*np.shape(heatmat)[1])+xextra, min(maxsize, 0.3*np.shape(heatmat)[0])+yextra), dpi = dpi)
    if 0.3*np.shape(heatmat)[1] > maxsize*ratio[0]/ratio[1]:
        xticklabels = None
        plot_value = False
    if 0.3*np.shape(heatmat)[0] > maxsize:
        yticklabels = None
        plot_value = False
    
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_position([0.15,0.15,0.7,0.75])
    ax.tick_params(which = 'both', bottom = False, labelbottom = False, left = False, labelleft = False)
    
    
    if sortx is not None:
        axdenx = fig.add_subplot(711)
        axdenx.spines['top'].set_visible(False)
        axdenx.spines['right'].set_visible(False)
        axdenx.spines['bottom'].set_visible(False)
        axdenx.tick_params(which = 'both', bottom = False, labelbottom = False)
        axdenx.set_position([0.15,0.91,0.7,0.19])
        dnx = dendrogram(Zx, ax = axdenx, no_labels = True, above_threshold_color = 'k', color_threshold = color_cutx, orientation = 'top')
        
        sortx = dnx['leaves']
        heatmat = heatmat[:, sortx]
        if x_attributes is not None:
            x_attributes = x_attributes[:, sortx]
            
        if xticklabels is not None:
            xticklabels = xticklabels[sortx]
            
        if xdenline is not None:
            axdenx.plot([0,len(heatmat[0])*10], [xdenline, xdenline], color = 'r')
    else:
        sortx = np.arange(len(heatmat[0]), dtype = int)
    
    sys.setrecursionlimit(100000)    
    
    if sorty =='maxdist':
        sorty = sortsize
    elif sorty is not None:
        axdeny = fig.add_subplot(171)
        axdeny.spines['top'].set_visible(False)
        axdeny.spines['right'].set_visible(False)
        axdeny.spines['left'].set_visible(False)
        axdeny.tick_params(which = 'both', left = False, labelleft = False)
        axdeny.set_position([0.05,0.15,0.09,0.75])
        dny = dendrogram(Zy, ax = axdeny, no_labels = True, color_threshold = color_cuty, above_threshold_color = 'k', orientation = 'left', get_leaves = True)
        sorty = dny['leaves']
        heatmat = heatmat[sorty]
        #axdeny.set_yticks(axdeny.get_yticks()[1:])

        if y_attributes is not None:
            y_attributes = y_attributes[sorty]
            
        if yticklabels is not None:
            yticklabels = yticklabels[sorty]
        if ydenline is not None:
            axdeny.plot([ydenline, ydenline], [0,len(heatmat)*10], color = 'r')
    else:
        sorty = np.arange(len(heatmat), dtype = int)
    
    
    if vmin is None:
        vmin = np.amin(heatmat)
    if vmax is None:
        vmax = np.amax(heatmat)
    
    ax.imshow(heatmat, aspect = 'auto', cmap = heatmapcolor, vmin = vmin, vmax = vmax, origin = 'lower')
    ax.set_yticks(np.arange(len(heatmat)))
    ax.set_xticks(np.arange(len(heatmat[0])))
    
    if plot_value:
        if np.amax(np.absolute(heatmat)) >= 10:
            heattext = np.array(heatmat, dtype = int)
        else:
            heattext = np.around(heatmat, 2)
        for c in range(len(heattext[0])):
            for d in range(len(heattext)):
                ax.text(c,d,str(heattext[d,c]), color = 'k', ha = 'center', fontsize = 6)
    
    
    if grid:
        ax.set_yticks(np.arange(len(heatmat)+1)-0.5, minor = True)
        ax.set_xticks(np.arange(len(heatmat[0])+1)-0.5, minor = True)
        ax.grid(color = 'k', which = 'minor')


    if x_attributes is not None:
        for x, xunique in enumerate(x_attributes):
            xunique = np.unique(xunique)
            for s, xuni in enumerate(xunique):
                x_attributes[x, x_attributes[x] == xuni] = s
        x_attributes = x_attributes.astype(int)
        axatx = fig.add_subplot(717)
        axatx.spines['top'].set_visible(False)
        axatx.spines['bottom'].set_visible(False)
        axatx.spines['right'].set_visible(False)
        axatx.spines['left'].set_visible(False)
        axatx.tick_params(which = 'both', bottom = False, labelbottom = False, left = False, labelleft = False, labelright = True)
        axatx.set_position([0.15,0.11,0.7,0.04])
        axatx.imshow(x_attributes, aspect = 'auto', cmap = xatt_color)
        axatx.set_xticks(np.arange(len(heatmat[0])))        
        if xlabel is not None:
            axatx.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False)
            axatx.set_xlabel(xlabel)
        if xticklabels is not None:
            axatx.tick_params(which = 'both', bottom = True, labelbottom = True, left = False, labelleft = False)
            axatx.set_xticklabels(xticklabels, rotation  = 90)
    
    elif xlabel is not None:
        ax.set_xlabel(xlabel)
    elif xticklabels is not None:
        ax.tick_params(which = 'both', bottom = True, labelbottom = True, left = False, labelleft = False)
        ax.set_xticklabels(xticklabels, rotation = 90)
        
    
    if y_attributes is not None:
        for y, yunique in enumerate(y_attributes.T):
            yunique = np.unique(yunique)
            for s, yuni in enumerate(yunique):
                y_attributes[y_attributes[:,y] == yuni,y] = s
        y_attributes = y_attributes.astype(int)
        axaty = fig.add_subplot(177)
        axaty.spines['top'].set_visible(False)
        axaty.spines['bottom'].set_visible(False)
        axaty.spines['right'].set_visible(False)
        axaty.spines['left'].set_visible(False)
        axaty.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False)
        axaty.set_position([0.85,0.15,0.03,0.75])
        axaty.imshow(y_attributes, aspect = 'auto', cmap = yatt_color)
        axaty.set_yticks(np.arange(len(heatmat)))
        if ylabel is not None:
            axaty.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False, labelright = True)
            axaty.set_ylabel(ylabel)
        if yticklabels is not None:
            axaty.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False, labelright = True, right = True)
            axaty.set_yticklabels(yticklabels)
    elif ylabel is not None:
        if sorty is not None:
            axdeny.set_ylabel(ylabel)
        else:
            ax.set_ylabel(ylabel)
    elif yticklabels is not None:
        if xticklabels is None:
            ax.tick_params(which = 'both', bottom = False, labelbottom = False, left = False, labelleft = False, labelright = True)
        else:
            ax.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False, labelright = True)
        ax.set_yticklabels(yticklabels)
    else:
        ax.set_yticks(np.arange(0, len(heatmat),200))
        if xticklabels is None:
            ax.tick_params(which = 'both', bottom = False, labelbottom = False, left = False, labelleft = False, labelright = True)
        else:
            ax.tick_params(which = 'both', bottom = False, labelbottom = True, left = False, labelleft = False, labelright = True)
        
    

    if figname is not None:
        fig.savefig(figname+'_heatmap.jpg', bbox_inches = 'tight')
        print( 'SAVED', figname)
    else:
        plt.show()
    plt.close()
    return sortx, sorty

# Function to plot the distribution of the number or SNVs per person
def plot_distribution(nsnps, xlabel, ylabel = 'Number', bins = 40, color = 'indigo', alpha = 0.8, histtype = 'bar', savefig = None, fmt = '.jpg'):
    fig = plt.figure(figsize = (4,3.5), dpi = 200)
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.hist(nsnps, bins = bins, color = color, alpha = alpha, histtype = histtype)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if savefig is not None:
        fig.savefig(os.path.splitext(savefig)[0]+'_distribution'+fmt, dpi = 300, bbox_inches = 'tight')
    else:
        fig.tight_layout()
        plt.show()




if __name__ == '__main__': 
    # Read SNP file
    snps = np.load(sys.argv[1])
    snames = snps['rows']
    inames = snps['cols']
    snps = snps['snps']
    
    # select set of individual that was used in analysis
    if '--select_individuals' in sys.argv:
        indsel = np.genfromtxt(sys.argv[sys.argv.index('--select_individuals')+1], dtype = str)
        mask = np.isin(inames, indsel)
        snps, inames = snps[:, mask], inames[mask]
    outname = None

    if '--savefig' in sys.argv:
        outname = sys.argv[sys.argv.index('--savefig')+1]
    print(np.shape(snps))
    plot_distribution(np.sum(snps>0, axis = 0), 'Number SNVs', ylabel = 'Number individuals', savefig = outname)
    
    # only look at snv presence not at individual alleles
    if '--snv_presence' in sys.argv:
        snps = (snps > 0).astype(float)
        print(snps)


    # filter snps with low frequency
    snpfreq = np.sum(snps > 0, axis = 1)/len(snps[0])
    if '--minfreq' in sys.argv:
        minf = float(sys.argv[sys.argv.index('--minfreq')+1])
        mask = snpfreq > minf
        snames = snames[mask]
        snps = snps[mask]
        print('minfreq', minf, np.shape(snps))
    
    # combine snps that correleate more than x with each other (also referred to as in LD)
    if '--combine_snps' in sys.argv:
        dthrsh = 1.-float(sys.argv[sys.argv.index('--combine_snps')+1])
        clu = agc(n_clusters = None, affinity = 'correlation', linkage = 'complete', distance_threshold = dthrsh).fit(snps)
        labels = clu.labels_
        nsnps, nsnames = [], []
        for u, ul in enumerate(np.unique(labels)):
            mask = labels == ul
            nsnps.append(np.mean(snps[mask],axis = 0))
            nname = '-'.join(snames[mask])
            if len(nname) > 50:
                nname = 'C('+snames[mask][0]+')'
            nsnames.append(nname)
        snps, snames = np.array(nsnps), np.array(nsnames)
        print('combine snps', dthrsh, np.shape(snps))
    
    # minimum number of snps per individual
    if '--minsnps' in sys.argv:
        msnp = int(sys.argv[sys.argv.index('--minsnps')+1])
        mask = np.sum(snps,axis = 0) > msnp
        print(np.sum(mask))
        snps = snps[:,mask]
        inames = inames[mask]
        print('minsnps', msnp, np.shape(snps))
    
    # combine individuals if they have very similar genotype
    if '--combine_individuals' in sys.argv:
        dthrsh = 1.-float(sys.argv[sys.argv.index('--combine_individuals')+1])
        clu = agc(n_clusters = None, affinity = 'correlation', linkage = 'complete', distance_threshold = dthrsh).fit(snps.T)
        labels = clu.labels_
        nsnps, ninames = [], []
        for u, ul in enumerate(np.unique(labels)):
            mask = labels == ul
            nsnps.append(np.mean(snps[:,mask],axis = 1))
            nname = '-'.join(inames[mask])
            if len(nname) > 50:
                nname = 'C('+inames[mask][0]+')'
            ninames.append(nname)
        snps, inames = np.array(nsnps).T, np.array(ninames)
        print('combine individuals', dthrsh, np.shape(snps))


    plot_heatmap(snps, measurex = 'euclidean', measurey = 'euclidean', sortx = 'average', sorty = 'average', x_attributes = None, y_attributes = None, xattr_name = None, yattr_name = None, heatmapcolor = cm.Purples, xatt_color = None, yatt_color = None, pwms = None, combine_cutx = 0., combine_cuty = 0., color_cutx = 0., color_cuty = 0., plot_value = False, vmin = None, vmax = None, grid = True, xdenline = None, ydenline = None, xlabel = 'Individuals', ylabel = None, xticklabels = inames, yticklabels  = snames, dpi = 100, figname = outname, maxsize = 30)

    

 




