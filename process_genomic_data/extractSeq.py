import pandas as pd
import numpy as np
import gzip
import sparse
import sys

pittPath = '/bgfs/mchikina/byng/'
filepath = pittPath+'rosmapAD/projects/insilicoMutagenesis/extractSequence/results/'
refpath = pittPath+'humanRefGenome/data/hg38/'

# Parameters
#ch = 22 # Select chromosome
ch = int(sys.argv[1])
gList = sys.argv[2]
nSym = 8 # ATGCx2 for paternal and maternal
nSub = 1161
win = 1e5
if win==1e4:
    winSiz = '10K'
elif win==1e5:
    winSiz = '100K'
wpl = 50 # Number of words per line in hg38

# Load gene names and windows
#data = pd.read_csv(filepath+'geneWin'+winSiz+'.txt',sep='\t',header=None)
data = pd.read_csv(filepath+'geneWin'+winSiz+gList+'.txt',sep='\s+',header=None)
data.columns = ['ensg','chr','winS','winE']
#data = data.sort_values(by='chr',axis=0)
data = data.loc[data.chr==ch]
nGene = len(data.index)    
print('geneWin'+winSiz+gList) # For checking if all genes ran

# Loop through genes
for j in np.arange(nGene):
    print(j)   
    
    # Load reference sequence
    with gzip.open(refpath+'chr'+str(ch)+'.fa.gz','rt') as f: 
        f.readline() # Remove non-sequence line
        #ref = f.read().replace('\n','') # Read file as a single string with \n removed but uses too much RAM
        start = np.mod(data.iloc[j,2],wpl)-1
        nLine = np.floor_divide(data.iloc[j,2],wpl)
        if start==-1: # gene j at the end of a line
            nLine -= 1  
        for l in np.arange(nLine):
            f.readline() # Get to line where gene j is located
        if start==-1:
            f.read(wpl-1) # Get to end of the current line
            ref = f.read(int(2*win+1+np.floor_divide(2*win+1,wpl)+1)) # Read 2*win+1 bases + \n's + 1
        else:
            f.read(start) # Get to start of gene j location
            ref = f.read(int(2*win+1+np.floor_divide(2*win+1,wpl))) # Read 2*win+1 bases + \n's
    
    # Extract reference sequence within window
    #seqR = np.array(list(ref[data.iloc[j,2]-1:data.iloc[j,3]].upper()))
    seqR = np.array(list(ref.replace('\n','').upper()))

    # Check if variants exist
    f = open(filepath+'variantNucleotide'+winSiz+'/'+data.iloc[j,0]+'.csv')
    line = f.readline()
    f.close()
    if line != '': 
        # Load variant nucleotide
        snp = pd.read_csv(filepath+'variantNucleotide'+winSiz+'/'+data.iloc[j,0]+'.csv',header=None,encoding='latin1')

        # Extract paternal variants
        snpP = snp.replace(['A|A','A|T','A|G','A|C'],'A')
        snpP = snpP.replace(['T|A','T|T','T|G','T|C'],'T')
        snpP = snpP.replace(['G|A','G|T','G|G','G|C'],'G')
        snpP = snpP.replace(['C|A','C|T','C|G','C|C'],'C')

        # Extract maternal variants
        snpM = snp.replace(['A|A','T|A','G|A','C|A'],'A')
        snpM = snpM.replace(['A|T','T|T','G|T','C|T'],'T')
        snpM = snpM.replace(['A|G','T|G','G|G','C|G'],'G')
        snpM = snpM.replace(['A|C','T|C','G|C','C|C'],'C')

        # Insert variants to reference sequence
        nVar = len(snp.index)

        # Loop over subjects
        onehot = np.zeros((nSym,int(2*win+1),nSub),dtype='i8')
        for k in np.arange(nSub):
            seqP = seqR.copy()
            seqM = seqR.copy()

            # Loop over variants
            for i in np.arange(nVar): 
                seqP[snpP.iloc[i,1]-data.iloc[j,2]] = snpP.iloc[i,k+2]
                seqM[snpM.iloc[i,1]-data.iloc[j,2]] = snpM.iloc[i,k+2]

            # Convert to one-hot encoding
            onehot[:,:,k] = [seqP=='A',seqP=='T',seqP=='G',seqP=='C',seqM=='A',seqM=='T',seqM=='G',seqM=='C']
    else:
        print('no variants')
        
        # Loop over subjects
        onehot = np.zeros((nSym,int(2*win+1),nSub),dtype='i8')
        for k in np.arange(nSub):        
            # Convert to one-hot encoding
            onehot[:,:,k] = [seqR=='A',seqR=='T',seqR=='G',seqR=='C',seqR=='A',seqR=='T',seqR=='G',seqR=='C']
        
    # Convert to sparse and save
    onehot = sparse.COO.from_numpy(onehot)
    sparse.save_npz(filepath+'/sequence'+winSiz+'/chr'+str(ch)+'/'+data.iloc[j,0],onehot)
