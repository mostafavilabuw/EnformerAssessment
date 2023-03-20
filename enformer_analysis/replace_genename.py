# Replace gene names with gene IDs

import numpy as np
import sys, os

a = np.genfromtxt('Full_analysis/gene-ids-and-positions.tsv', dtype = str, delimiter = '\t', skip_header = 1)

repfile = open(sys.argv[1],'r').readlines()
obj = open(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'ENSG.tsv', 'w')
names = list(a[:,0])
for h, head in enumerate(repfile):
    oldname = head.split()[0]
    if oldname in names or h > :
        newname = a[names.index(oldname),1]
        print(oldname, newname)
        obj.write(head.replace(oldname, newname))
    else:
        obj.write(head)

