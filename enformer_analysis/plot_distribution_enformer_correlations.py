# Plot the distribution of genes with a SuSie assigned SNP versus the ones without
#python3 plot_distribution_enformer_correlations.py Full_analysis/Prediction_correlationsCageAdultBrain.txt susie_SNP_gene_CortexENSG.txt

import numpy as np
import sys, os
import glob
import matplotlib.pyplot as plt
import seaborn as sns

def vplot(corrs, ylabel = None, xlabel = None, names = None, outname = None):
    medians = []
    quartile1 = []
    quartile3 = []
    for f, corr in enumerate(corrs):
        medians.append(np.median(corr))
        quartile1.append(np.percentile(corr, 25))
        quartile3.append(np.percentile(corr, 75))
    fig = plt.figure(figsize = (len(corrs)*0.8,3.5))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.violinplot(data = corrs, ax = ax, width = 1., alphas = 0.5, palette = 'magma', colors = [0.1, 1.2], vmin = 0, vmax = 2, cut = 0)
    #sns.swarmplot(data = corrs, ax = ax, color = 'k', size = 3, zorder = 2)
    #ax.vlines(np.arange(len(corrs)), quartile1, quartile3, color='silver', linestyle='-', lw=1, zorder = 3)
    #ax.scatter(np.arange(len(corrs)), medians, color = 'silver',zorder = 3)
    ax.set_xticks(np.arange(len(corrs)))
    if names is not None:
        ax.set_xticklabels(names)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if outname is not None:
        fig.savefig(outname, transparent = True, dpi = 500, bbox_inches = 'tight')
        print(outname)
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':
    # Read correlation to observed data
    corr = np.genfromtxt(sys.argv[1], dtype = str)
    cnames, corr = corr[:,0], corr[:,1].astype(float)
    cnames, corr = cnames[~np.isnan(corr)], corr[~np.isnan(corr)]
    # Read list of gene set
    clist = np.genfromtxt(sys.argv[2], dtype = str)[:,0]

    corrs = [corr[np.isin(cnames, clist)], corr[~np.isin(cnames, clist)]]
    print('neg', np.sum(corrs[0]<0))



    vplot(corrs, ylabel = 'R Enformer to obs.', names = ['Susie','Non-Susie'], outname = 'Distribution_EnformerCorrelation_susieset.jpg')
    vplot([np.absolute(corr) for corr in corrs], ylabel = 'Abs(R) Enformer to obs.', names = ['Susie','Non-Susie'], outname = 'Distribution_EnformerAbsCorrelation_susieset.jpg')


