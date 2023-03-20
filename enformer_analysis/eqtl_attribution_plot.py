# plot attributions versus eqtl values

import numpy as np
import sys, os
import matplotlib.pyplot as plt

# Read attribution values for all SNVs
attributions = np.genfromtxt(sys.argv[1])
# Read eqtls
eqtl = np.nan_to_num(np.genfromtxt(sys.argv[2]))
# sort if necessary
sort = np.argsort(attributions[:,0])[np.isin(np.sort(attributions[:,0]), eqtl[:,0])]
attributions = attributions[sort]
sort = np.argsort(eqtl[:,0])[np.isin(np.sort(eqtl[:,0]), attributions[:,0])]
eqtl = eqtl[sort]
#check sorting
if not np.array_equal(attributions[:,0], eqtl[:,0]):
    print('eqtl and attributions dont match')
    sys.exit()

xlabel = sys.argv[3]
ylabel = sys.argv[4]

if '--colors' in sys.argv:
    colors = np.genfromtxt(sys.argv[sys.argv.index('--colors')+1]) # snp_info file
    sort = np.argsort(colors[:,0])[np.isin(np.sort(colors[:,0]), attributions[:,0])]
    colors = colors[sort]
    if not np.array_equal(attributions[:,0], colors[:,0]):
        print('colors and attributions dont match')
        sys.exit()
    colors = colors[:,1]
    # select SNVs based on color assignment
    if '--minfreq' in sys.argv:
        mask = colors > float(sys.argv[sys.argv.index('--minfreq')+1])
        attributions = attributions[mask]
        eqtl = eqtl[mask]
        colors = colors[mask]

stdatt = np.std(attributions[:,1])
stdeqtl = np.std(eqtl[:,1])

maxatt = np.amax(np.absolute(attributions[:,1]))
if '--norm_attribution' in sys.argv:
    attributions[:,1] /= maxatt
    stdatt /= maxatt
    maxatt = 1
print(stdatt, stdeqtl)

fig = plt.figure(figsize = (4.,3.5))
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.plot([0,0],[-1,1], c = 'silver', lw = 0.7)
ax.plot([-maxatt,maxatt],[0,0], c = 'silver', lw = 0.7)
#ax.plot([stdatt,stdatt],[-1,1], c = 'green', lw = 0.7, ls = '--')
#ax.plot([-stdatt,-stdatt],[-1,1], c = 'green', lw = 0.7, ls = '--')
#ax.plot([-maxatt,maxatt],[stdeqtl,stdeqtl], c = 'green', lw = 0.7, ls = '--')
#ax.plot([-maxatt,maxatt],[-stdeqtl,-stdeqtl], c = 'green', lw = 0.7, ls = '--')
sort = np.argsort(colors)
cab = ax.scatter(attributions[sort,1], eqtl[sort,1], cmap = 'Blues', vmin = 0, vmax =1, c = colors[sort], edgecolor = 'grey')
fig.colorbar(cab, pad = 0., fraction = 0.09, shrink = 0.15, aspect = 2, anchor = (0.,0.99))

# Mark driver SNVs in plot with enlarged dots and red edgecolors
if '--drivers' in sys.argv:
    dobj = open(sys.argv[sys.argv.index('--drivers')+1], 'r').readlines()
    if len(dobj) > 0:
        dfile = [line.strip().split() for line in dobj]
        dfile = np.array(dfile,dtype =float)
        dfile = dfile[np.isin(dfile[:,0], attributions[:,0])]
        dloc = np.where(np.isin(attributions[:,0], dfile[:,0]))
        dloc = dloc[0]
        size0 = plt.rcParams['lines.markersize'] ** 2
        ax.scatter(attributions[dloc,1], eqtl[dloc,1], s = size0 * (1+2*dfile[:,-3]), linewidths = 1.2, cmap = 'Blues', vmin = 0, vmax =1, c = colors[dloc], edgecolor = 'red')
        dtypefile = open(os.path.splitext(sys.argv[sys.argv.index('--drivers')+1])[0]+'_types.txt', 'w')
        for t in dloc:
            print(int(attributions[t,0]), attributions[t,1])
            if '--name_driver' in sys.argv:
                ax.text(int(attributions[t,1]), eqtl[t,1], str(int(attributions[t,0])), ha = 'left')
            # Determine SNP type of drivers by looking at position in eqtl and attribution plot
            # For each driver determine where it's located based on the std of eqtls and attributions
            dtypefile.write(str(int(attributions[t,0]))+' '+str(round(attributions[t,1]/stdatt,2))+' '+str(round(eqtl[t,1]/stdeqtl,2))+'\n')

ax.set_xticks([-np.around(maxatt,3), 0, np.around(maxatt,3)])

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
if '--show' in sys.argv:
    plt.show()
else:
    dpi, fmt = 250, '.jpg'
    if '--dpi' in sys.argv:
        dpi = int(sys.argv[sys.argv.index('--dpi')+1])
    if '--fmt' in sys.argv:
        fmt = sys.argv[sys.argv.index('--fmt')+1]
    fig.savefig(os.path.splitext(sys.argv[1])[0]+'_vs_'+sys.argv[5]+fmt, transparent = True, dpi = dpi, bbox_inches = 'tight')




