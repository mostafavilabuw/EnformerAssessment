# Generate file that computes the distance of the drivers to the TSS

import numpy as np
import sys, os
import glob 

dfiles = np.sort(glob.glob('ENSG*'+sys.argv[1]))
tssfile = np.genfromtxt(sys.argv[2], dtype = str)

for d, df in enumerate(dfiles):
    gene = df.strip(sys.argv[1])
    dobj = open(df, 'r').readlines()
    if len(dobj) > 0:
        dfile = [line.strip().split() for line in dobj]
        dfile = np.array(dfile,dtype =float)
        tss = int(tssfile[list(tssfile[:,0]).index(gene),1])
        for e, var in enumerate(dfile):
            print(gene, int(var[0]), int(var[0]-tss), var[-3])



