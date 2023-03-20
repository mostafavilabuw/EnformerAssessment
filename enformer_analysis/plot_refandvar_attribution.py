# Plot the mean attribution and attribution for base to base changes in a specific window
#python3 ../plot_refandvar_attribution.py ENSG00000134202_109741038_attribution.npz --squaresize 0.12 --include_snvs 109741038 ../ism_res/ENSG00000134202_ism_attributions.txt ../eqtl/ENSG00000134202_corr.txt --markdrivers ../ism_res/ENSG00000134202_ism_attributions_driversfw.txt --dpi 350 --include_conservation ../PhyloP100/ENSG00000134202_1000tss_in_hg38.phyloP100way.txt --excludeheatmap --enlargepwm 1.8

import numpy as np
import matplotlib.pyplot as plt
import logomaker as lm
import sys, os
import pandas as pd

nts = list('ACGT')
ntsar = np.array(nts)

# Load npz with attributions
f = np.load(sys.argv[1])
outname = os.path.splitext(sys.argv[1])[0]
gene, locge = os.path.split(outname)[1].split('_')[:2]
attribution, seq = f['attribution'], str(f['seq'])
ohseq = np.array(list(seq))[:,None] == ntsar

# Norm attribuitons with external values, f.e. the max or a std
if '--norm_attributions' in sys.argv:
    norms = np.genfromtxt(sys.argv[sys.argv.index('--norm_attributions')+1], dtype = str)
    ngene = list(norms[:,0]).index(sys.argv[sys.argv.index('--norm_attributions')+2])
    norm = float(norms[ngene,int(sys.argv[sys.argv.index('--norm_attributions')+3])])
    attribution/=norm
maxg = np.amax(np.absolute(attribution))

# adjust figsize
if '--figsize' in sys.argv:
    l,h = sys.argv[sys.argv.index('--figsize')+1].split(',')
    fig = plt.figure(figsize = (float(l), float(h)))
    sqsize = float(l)/len(attribution)
# or adjust size per unit of attribution 
elif '--squaresize' in sys.argv:
    sqsize = float(sys.argv[sys.argv.index('--squaresize')+1])
    fig = plt.figure(figsize = (sqsize*len(attribution),sqsize*9))
else:
    fig = plt.figure(figsize = (0.4*len(attribution),3.5))
    sqsize = 0.4

# Plot the attributions of individual SNVs into the mean attribution plot
if '--include_snvs' in sys.argv:
    center = int(sys.argv[sys.argv.index('--include_snvs')+1])
    snvs = np.genfromtxt(sys.argv[sys.argv.index('--include_snvs')+2], dtype = str)
    freqs = np.genfromtxt(sys.argv[sys.argv.index('--include_snvs')+3], dtype = str)
    loc = snvs[:,0].astype(int)-center + int(len(attribution)/2)
    lmask = (loc < len(attribution))*(loc >=0 )
    loc, snvatt, snvs = loc[lmask], snvs[lmask,1].astype(float), snvs[lmask]
    fmask = (freqs[:,0].astype(int)-center + int(len(attribution)/2) < len(attribution)) * (freqs[:,0].astype(int)-center + int(len(attribution)/2) >= 0)
    if not np.array_equal(freqs[fmask, 0], snvs[:,0]):
        print('Freqs and attributions dont match')
        print(freqs[fmask, 0], snvs[:,0])
        sys.exit()
    freqs = np.nan_to_num(freqs[fmask,1].astype(float))
    if '--norm_attributions' in sys.argv:
        snvatt /=norm 
    # also want frequency and common variant to show in the heatmap
    edgecolors = np.chararray(len(snvatt), itemsize = 10, unicode = True)
    edgecolors[:] = 'k'
    msizes = np.ones(len(snvatt)) * sqsize * 450
    if '--markdrivers' in sys.argv:
        driverfile = open(sys.argv[sys.argv.index('--markdrivers')+1], 'r').readlines()
        if len(driverfile) > 0:
            drivers = np.array([line.strip().split() for line in driverfile])
            maindriver = np.argmax(drivers[:,-3].astype(float))
            driversinside = np.where(np.isin(drivers[:,0], snvs[:,0]))[0]
            if len(driversinside) > 0:
                edgecolors[np.isin(snvs[:,0], drivers[:,0])] = 'tomato'
                msizes[np.isin(snvs[:,0], drivers[:,0])] *= 2.5
                if maindriver in driversinside:
                    edgecolors[np.isin(snvs[:,0], drivers[maindriver,0])] = 'r'
# Include a plot with conservation values 
naxes = 2
if '--include_conservation' in sys.argv:
    consfile = np.genfromtxt(sys.argv[sys.argv.index('--include_conservation')+1], dtype = str)
    clocs, cons = consfile[:,1].astype(int)-center+ int(len(attribution)/2), consfile[:,2].astype(float)
    clmask = (clocs>=0)*(clocs<len(attribution))
    clocs, cons = clocs[clmask], cons[clmask]
    naxes = 3
    if '--excludeheatmap' in sys.argv:
        naxes = 2
    axcons = fig.add_subplot(naxes,1,2)
    axcons.spines['top'].set_visible(False)
    axcons.spines['right'].set_visible(False)
    axcons.spines['bottom'].set_visible(False)
    colors = np.chararray(len(cons), itemsize = 15, unicode = True)
    colors[:] = 'grey'
    #colors[cons<-1.3] = 'red'
    #colors[cons>1.3] = 'blue'
    cons[cons < -0.5] = -0.5
    conspwm = np.zeros(np.shape(attribution))
    conspwm[ohseq] = cons
    conspwm[conspwm<1.3] = 0
    cons[cons>=1.3] = 0
    axcons.bar(clocs, cons, color = colors)
    conspwm = pd.DataFrame({'A':conspwm[:,0],'C':conspwm[:,1], 'G':conspwm[:,2], 'T':conspwm[:,3]})
    lm.Logo(conspwm, ax = axcons)
    axcons.plot([0,len(attribution)], [0,0], color = 'grey')
    axcons.set_xlim([-0.5,len(attribution)-0.5])
    axcons.set_ylabel('PhyloP')

# Exclude the heatmap that shows attribuitons for individual base to base changes
if not '--excludeheatmap' in sys.argv:
    ax = fig.add_subplot(naxes,1,naxes)
    hm = ax.imshow(attribution.T, aspect = 'auto', vmin = -maxg, vmax = maxg, cmap = 'RdBu_r')
    ax.set_xticks(np.arange(len(attribution)))
    ax.set_yticks(np.arange(len(nts)))
    ax.set_yticklabels(nts)
    ax.set_xticklabels(list(seq))
    ax.set_xticks(np.arange(len(attribution))+0.5, minor = True)
    ax.set_yticks(np.arange(len(nts))+0.5, minor = True)
    ax.grid(which = 'minor', color = 'k')
    #fig.colorbar(hm, pad = 0., fraction = 0.09, shrink = 0.15, aspect = 2, anchor = (0.,0.99), ax = ax)
    ax.set_ylabel('ISM')


axpwm =  fig.add_subplot(naxes, 1, 1)
# Increase the height of the mean attribuition that is shown as a pwm
if '--enlargepwm' in sys.argv:
    pos1 = axpwm.get_position()
    enl = float(sys.argv[sys.argv.index('--enlargepwm')+1])
    axpwm.set_position([pos1.x0,pos1.y0,pos1.width,enl])
axpwm.set_title(gene+' - '+locge)
axpwm.spines['top'].set_visible(False)
axpwm.spines['right'].set_visible(False)
axpwm.spines['bottom'].set_visible(False)
axpwm.tick_params(bottom = False, labelbottom = False)
pwm = np.zeros(np.shape(attribution))
pwm[ohseq] = -np.mean(attribution, axis = 1)
axpwm.set_ylabel('Mean Attribution')
pwm = pd.DataFrame({'A':pwm[:,0],'C':pwm[:,1], 'G':pwm[:,2], 'T':pwm[:,3]})
lm.Logo(pwm, ax = axpwm)
ylim = axpwm.get_ylim()
if '--include_snvs' in sys.argv:
    axpwm.vlines(loc, np.zeros(len(snvatt)), snvatt, color = 'k', zorder = -1)
    axpwm.scatter(loc, snvatt, cmap = 'bwr', c = freqs, vmin = -0.5, vmax =0.5, marker = 'o', s = msizes, edgecolor = edgecolors, lw = 1.5)
    axpwm.set_ylim([min(ylim[0],np.amin(snvatt))*1.08,max(ylim[1],np.amax(snvatt))*1.08])

axpwm.plot([int(len(attribution)/2), int(len(attribution)/2)],ylim, c= 'darkgoldenrod', ls = '--')
if '--setylimpwm' in sys.argv:
    ylim = sys.argv[sys.argv.index('--setylimpwm')+1].split(',')
    axpwm.set_ylim([float(ylim[0]), float(ylim[1])])

dpi = 150
if '--dpi' in sys.argv:
    dpi = float(sys.argv[sys.argv.index('--dpi')+1])

fig.savefig(outname+'.jpg', dpi = dpi, bbox_inches = 'tight')
print(outname+'.jpg')
#plt.show()



