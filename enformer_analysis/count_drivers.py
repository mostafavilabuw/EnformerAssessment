import numpy as np
import sys, os 
import glob

ending = sys.argv[1]
files = np.sort(glob.glob('ENSG*'+ending))

for f, fi in enumerate(files):
    flen = open(fi,'r').readlines()
    print(fi.strip(ending), len(flen))

