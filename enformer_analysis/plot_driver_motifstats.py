# plots statistics of attributions around driver snps in comparison to global ism attributions
#python3 ../plot_driver_motifstats.py ../ism_res/ALLgenes_ism_attributions_driversfw_types.txt ALLgenes_ism_attributions_driversfw_types_ism_significance_stats.txt 2,15 --absolute --combine_columns max --lim 10 --cumulative -1 --nbins 22 --savefig ISMstats_maxzscore_in_var_or_ref.jpg

import numpy as np
import sys, os
import matplotlib.pyplot as plt

# File with drivers for all genes
drivers = np.genfromtxt(sys.argv[1], dtype = str)
# File generated with extract_ism_stats_around_drivers.py
motifs = np.genfromtxt(sys.argv[2], dtype = str)
drnames = np.array([m[0]+m[1] for m in drivers])
motnames = np.array([m[0]+m[1] for m in motifs])
# sort files to match SNVs
dsort = np.argsort(drnames)[np.isin(np.sort(drnames), motnames)]
drnames, drivers = drnames[dsort], drivers[dsort]
msort = np.argsort(motnames)[np.isin(np.sort(motnames), drnames)]
motnames, motifs = motnames[msort], motifs[msort]

# If selected, only main drivers will be selected
if '--maindrivers' in sys.argv:
    mainfile = np.genfromtxt(sys.argv[sys.argv.index('--maindrivers')+1], dtype = str)
    mcol = int(sys.argv[sys.argv.index('--maindrivers')+2])
    mains = []
    for g, gen in enumerate(np.unique(mainfile[:,0])):
        mask = np.where(mainfile[:,0] == gen)[0]
        if len(mask) > 1:
            argm = np.argmax(mainfile[mask,mcol].astype(float))
            mains.append(mainfile[mask[argm],0]+mainfile[mask[argm],1])
        else:
            mains.append(mainfile[mask[0],0]+mainfile[mask[0],1])
    mask = np.isin(drnames, mains)
    print('From', len(drnames))
    drnames, drivers, motnames, motifs = drnames[mask], drivers[mask], motnames[mask], motifs[mask]
    print('to', len(drnames))
    #print(drnames)

# check if rows are aligned
if np.array_equal(drnames, motnames):
    
    # determine driver type
    driverdir = drivers[:,-2].astype(float)
    drloc = drivers[:,1].astype(int)
    drivertype = []
    for d, dr in enumerate(drivers):
        drivertype.append(np.sign(float(dr[-1])*float(dr[-2])))
    drivers = drivers[:,0]
    drivertype = np.array(drivertype)

    # Potentially combine the entries of two columns, take the mean, max, min of these columns
    distarg = None
    if '--combine_columns' in sys.argv:
        column = np.array(sys.argv[3].split(','), dtype = int)
        ctype = sys.argv[sys.argv.index('--combine_columns')+1]
        dist = motifs[:, column].astype(float)
        if '--absolute' in sys.argv:
            dist = np.absolute(dist)
        if ctype == 'min':
            # Decide if you would like to split the histogram between data points taken from the ref and var sequence
            if '--splitrefvar' in sys.argv:
                distarg = np.argmin(dist, axis =1)
            dist = np.amin(dist, axis =1)
        elif ctype == 'max':
            if '--splitrefvar' in sys.argv:
                distarg = np.argmax(dist, axis =1)
            dist = np.amax(dist, axis =1)
        elif ctype == 'mean':
            dist = np.mean(dist, axis =1)

    else:
        column = int(sys.argv[3])
        dist = motifs[:, column].astype(float)
    
    # transform the data more
    if '--absolute' in sys.argv:
        dist = np.absolute(dist)
    
    if '--scale' in sys.argv:
        dist *= float(sys.argv[sys.argv.index('--scale')+1])

    if '--lim' in sys.argv:
        lim = int(sys.argv[sys.argv.index('--lim')+1])
        dist[np.absolute(dist)>lim] = np.sign(dist[np.absolute(dist)>lim])*lim

    if '--log' in sys.argv:
        dist = np.sign(dist)*np.log10(1+np.absolute(dist))
    if '--symlog' in sys.argv:
        dist[np.absolute(dist)>1] = np.sign(dist[np.absolute(dist)>1])*np.log10(10+np.absolute(dist[np.absolute(dist)>1]))
    
    # define umber of bins
    nbins = 31
    if '--nbins' in sys.argv:
        nbins = int(sys.argv[sys.argv.index('--nbins')+1])
    rnge = np.amax(dist)-np.amin(dist)
    over = rnge/nbins/2
    if '--nooverhang' in sys.argv:
        bins = np.linspace(np.amin(dist), np.amax(dist),nbins)
    # determine limits of x axis
    elif '--xlim' in sys.argv:
        blim0, blim1 = sys.argv[sys.argv.index('--xlim')+1].split(',')
        bins = np.linspace(float(blim0), float(blim1), nbins)
    else:
        bins = np.linspace(np.amin(dist)-over, np.amax(dist)+over,nbins)
    print(bins)
    
    # define if histogram is cumulative
    cum = False
    if '--cumulative' in sys.argv:
        cum = int(sys.argv[sys.argv.index('--cumulative')+1])
    
    # define if density or numbers are shown in histogram
    dens = True
    if '--numbers' in sys.argv:
        dens = False
    
    htype = 'bar'
    if '--histtype' in sys.argv:
        htype =sys.argv[sys.argv.index('--histtype')+1]
    
    # print some stats
    if '--printout' in sys.argv:
        for dt in [1,-1]:
            print('Drivertype', dt)
            for d in np.where(drivertype == dt)[0]:
                print(motifs[d,0], motifs[d,1], dist[d])

    fig = plt.figure(figsize = (3.5, 3.5), dpi = 200)
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if distarg is not None:
        nc, bi = np.histogram(dist[(drivertype == 1)*(distarg == 0)], bins = bins, density = dens)
        nc = np.cumsum(nc[::cum])[::cum]*np.sum((drivertype == 1)*(distarg == 0))/np.sum(drivertype == 1)
        xc = np.mean(np.array([bi[:-1], bi[1:]]),axis = 0)
        ax.bar(xc, nc, alpha = 0.6, color = 'navy', label = 'Supportive ref',zorder = 0,width= 1)
        ncv, bi = np.histogram(dist[(drivertype == 1)*(distarg == 1)], bins = bins, density = dens)
        ncv = np.cumsum(ncv[::cum])[::cum]*np.sum((drivertype == 1)*(distarg == 1))/np.sum(drivertype == 1)
        xc = np.mean(np.array([bi[:-1], bi[1:]]),axis = 0)
        ax.bar(xc, ncv, alpha = 0.6, bottom = nc , color = 'slateblue', label = 'Supportive var',zorder = -2, width = 1)
        nc, bi = np.histogram(dist[(drivertype == -1)*(distarg == 0)], bins = bins, density = dens)
        nc = np.cumsum(nc[::cum])[::cum]*np.sum((drivertype == -1)*(distarg == 0))/np.sum(drivertype == -1)
        xc = np.mean(np.array([bi[:-1], bi[1:]]),axis = 0)
        ax.bar(xc, nc, alpha = 0.6, color = 'darkgoldenrod', label = 'Unsupportive ref',zorder = 1, width= 1)
        ncv, bi = np.histogram(dist[(drivertype == -1)*(distarg == 1)], bins = bins, density = dens)
        ncv = np.cumsum(ncv[::cum])[::cum]*np.sum((drivertype == -1)*(distarg == 1))/np.sum(drivertype == -1)
        xc = np.mean(np.array([bi[:-1], bi[1:]]),axis = 0)
        ax.bar(xc, ncv, alpha = 0.6, bottom = nc , color = 'goldenrod', label = 'Unsupportive var',zorder = -1, width = 1)
        
    else:
        nc, bi, t_ = ax.hist(dist[(drivertype == 1)], bins = bins, alpha = 0.6, color = 'navy', histtype = htype, density = dens, cumulative = cum, label = 'Supportive')
        print(nc)
        print(dist[(drivertype == 1)])
        nc, bi, t_ = ax.hist(dist[(drivertype == -1)], bins = bins, alpha = 0.6, color = 'goldenrod', histtype =htype, density = dens, cumulative = cum, label = 'Unsupportive')
        print(nc)
        print(dist[(drivertype == -1)])
    
    if '--xlabel' in sys.argv:
        ax.set_xlabel(sys.argv[sys.argv.index('--xlabel')+1])
    if '--ylabel' in sys.argv:
        ax.set_ylabel(sys.argv[sys.argv.index('--ylabel')+1])
    
    legloc = 'best'
    if '--legendpos' in sys.argv:
        legloc = sys.argv[sys.argv.index('--legendpos')+1]

    ax.legend(loc=legloc)
    if '--savefig' in sys.argv:
        fig.savefig(sys.argv[sys.argv.index('--savefig')+1], dpi = 250, bbox_inches = 'tight')
        
    else:
        plt.show()

