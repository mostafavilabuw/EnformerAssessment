import numpy as np
import sys, os
import glob


files=np.sort(glob.glob('ENSG*'+sys.argv[1]))

for f, fil in enumerate(files):
    gene = fil.strip(sys.argv[1])
    dobj = open(fil, 'r').readlines()
    if len(dobj) > 0:
        for line in dobj:
            line = line.strip().split()
            print(gene, line[0], line[1], line[2])

