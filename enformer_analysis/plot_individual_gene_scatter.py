# Plots scatter plot between observed and predicted expression values across same set of individuals for a single gene
# $python plot_individual_gene_scatter.py Observed_gene_expression.txt Enformer_predictions.txt DDX11 --figsize 4 3

import numpy as np
import sys, os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def read(file):
    genes, indv, exp = [],[],[]
    for l, line in enumerate(open(file,'r').readlines()):
        if l == 0:
            genes = line.strip().split('\t')
        else:
            line = line.strip().split('\t')
            indv.append(line[0])
            exp.append(line[1:])
    return np.array(genes), np.array(indv), np.array(exp, dtype = float)

# Read observed expression values
obgenes, obindv, obexp = read(sys.argv[1])
# Read predicted expression values
enfgenes, enfindv, enfexp = read(sys.argv[2])
# Define gene (column)
chg = sys.argv[3]

# Sort enfgenes and obgenes and the data matrices to match
e_, s1a = np.unique(enfgenes, return_index = True)
o_, s2a = np.unique(obgenes, return_index = True)
s1a, s1b = s1a[np.isin(e_, obgenes)], np.argsort(enfindv)[np.isin(np.sort(enfindv), obindv)]
s2a, s2b = s2a[np.isin(o_, enfgenes)], np.argsort(obindv)[np.isin(np.sort(obindv), enfindv)]
enfgenes, enfindv, enfexp = enfgenes[s1a], enfindv[s1b], enfexp[s1b][:,s1a]
obgenes, obindv, obexp = obgenes[s2a], obindv[s2b], obexp[s2b][:,s2a]

# Check if rows and columns match from sorting
print(np.array_equal(enfgenes, obgenes), np.array_equal(enfindv, obindv), np.shape(enfexp), np.shape(obexp))

figsize = (4,3.8)
if '--figsize' in sys.argv:
    figsize = (float(sys.argv[sys.argv.index('--figsize')+1]), float(sys.argv[sys.argv.index('--figsize')+2]))


# Select column of gene
ic = list(enfgenes).index(chg)
fig = plt.figure(figsize=figsize, dpi = 200)
ax = fig.add_subplot(111)
r = round(pearsonr(obexp[:,ic], enfexp[:,ic])[0], 2)
ax.scatter(obexp[:,ic], enfexp[:,ic], alpha = 0.8, label = 'R='+str(r),color = 'darkred')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Observed')
ax.set_ylabel('Predicted')
ax.legend()
ax.set_title(chg)
if '--setxlim' in sys.argv:
    ax.set_xlim(np.array(sys.argv[sys.argv.index('--setxlim')+1].split(','), dtype = float))

if '--setylim' in sys.argv:
    ax.set_ylim(np.array(sys.argv[sys.argv.index('--setylim')+1].split(','), dtype = float))

fig.savefig(chg+'scatter.jpg', dpi = 350, bbox_inches = 'tight')
if '--show' in sys.argv:
    plt.show()


