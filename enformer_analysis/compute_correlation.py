import numpy as np
import sys, os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def read(file, delimiter = ' '):
    genes, indv, exp = [],[],[]
    for l, line in enumerate(open(file,'r').readlines()):
        if l == 0:
            genes = line.strip('#').strip().split(delimiter)
        else:
            line = line.strip().split(delimiter)
            indv.append(line[0])
            exp.append(line[1:])
    return np.array(genes), np.array(indv), np.array(exp, dtype = float)

delimiter = ' '
if '--delimiter' in sys.argv:
    delimiter = sys.argv[sys.argv.index('--delimiter')+1]

obgenes, obindv, obexp = read(sys.argv[1], delimiter = delimiter)
enfgenes, enfindv, enfexp = read(sys.argv[2], delimiter = delimiter)

e_, s1a = np.unique(enfgenes, return_index = True)
o_, s2a = np.unique(obgenes, return_index = True)

s1a, s1b = s1a[np.isin(e_, obgenes)], np.argsort(enfindv)[np.isin(np.sort(enfindv), obindv)] 
s2a, s2b = s2a[np.isin(o_, enfgenes)], np.argsort(obindv)[np.isin(np.sort(obindv), enfindv)]

enfgenes, enfindv, enfexp = enfgenes[s1a], enfindv[s1b], enfexp[s1b][:,s1a]
obgenes, obindv, obexp = obgenes[s2a], obindv[s2b], obexp[s2b][:,s2a]
print(np.array_equal(enfgenes, obgenes), len(enfgenes), len(obgenes))


if '--get_allstats' in sys.argv:
    print('# Gene PearsonR P-value MeanObs StdObs MeanEnf StdEnf')
    for g, gene in enumerate(enfgenes):
        pears = pearsonr(enfexp[:,g], obexp[:,g])
        print(gene, round(pears[0],3), round(pears[1],3), round(np.mean(obexp[:,g]),2), round(np.std(obexp[:,g]),3), round(np.mean(enfexp[:,g]),2), round(np.std(enfexp[:,g]),3)) #, np.unique(enfexp[:,g])[0],  len(np.unique(enfexp[:,g])), len(np.unique(obexp[:,g])))
else:
    for g, gene in enumerate(enfgenes):
        pears = pearsonr(enfexp[:,g], obexp[:,g])
        print(gene, round(pears[0],3))



