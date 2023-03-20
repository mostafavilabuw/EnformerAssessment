# For all the drivers generate a file that contains the attribution normed by the max absolute attribution and their eqtl value

import numpy as np
import sys, os

attributions = np.genfromtxt(sys.argv[1])
eqtl = np.nan_to_num(np.genfromtxt(sys.argv[2]))

stdatt = np.amax(np.abs(attributions[:,1]))
stdeqtl = 1. #np.std(eqtl[:,1])
maxeqtl = np.amax(np.abs(eqtl[:,1]))
# sort attribution and eqtl
attributions = attributions[np.argsort(attributions[:,0])[np.isin(np.sort(attributions[:,0]), eqtl[:,0])]]
eqtl = eqtl[np.argsort(eqtl[:,0])[np.isin(np.sort(eqtl[:,0]), attributions[:,0])]]
if not np.array_equal(attributions[:,0], eqtl[:,0]):
    print('Not sorted correctly, check files')
    sys.exit()
    

dobj = open(sys.argv[3], 'r').readlines()
if len(dobj) > 0:
    dfile = [line.strip().split() for line in dobj]
    dfile = np.array(dfile,dtype =float)
    dloc = np.where(np.isin(attributions[:,0], dfile[:,0]))
    dloc = dloc[0]
    dtypefile = open(os.path.splitext(sys.argv[3])[0]+'_types.txt', 'w')
    for t in dloc:
        dtypefile.write(str(int(attributions[t,0]))+' '+str(round(attributions[t,1]/stdatt,3))+' '+str(round(eqtl[t,1]/stdeqtl,6))+'\n')





