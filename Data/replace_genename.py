import numpy as np
import sys, os


a = np.genfromtxt('gene-ids-and-positions.tsv', dtype = str, delimiter = '\t', skip_header = 1)
repfile = open(sys.argv[1],'r').readlines()
header = repfile[0].strip().split('\t')
print(header)
names = list(a[:,0])
nid = []
for h, head in enumerate(header):
    newname = a[names.index(head),1]
    if newname in nid:
        print(newname)
        print(header[nid.index(newname)], head)
    nid.append(newname)

obj = open(os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'ENSG.tsv', 'w')
obj.write('\t'+'\t'.join(np.array(nid))+'\n')
for r in repfile[1:]:
    obj.write(r)


