# Determine SNVs that drive the non-linear predictions by comparing partial sums of attributions times the genotype to the full predictions

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

# Read predictions for indivduals 
enfgenes, enfindv, enfexp = read('../Enformer_predictions_CageAdultBrain.txt')

# Read computed attributions
avals = np.genfromtxt(sys.argv[1], dtype = str)
avars = avals[:,0]
avals = avals[:,1].astype(float)

# Read genotypes for all individuals
snpfile = np.load(sys.argv[2])
indv = snpfile['columns'].astype(str)
snp = snpfile['snps']
locs = snpfile['rows'].astype(str)

# sort rows of genotypes to ism values
lsort = np.argsort(locs)[np.isin(np.sort(locs), avars)]
snp, locs = snp[lsort], locs[lsort]
asort = np.argsort(avars)[np.isin(np.sort(avars), locs)]
avals, avars = avals[asort], avars[asort]
if not np.array_equal(avars, locs):
    print('ATTENTION locs and avars not identical', len(avars), len(locs))
    sys.exit()

# Select the column of predicted expression values
gene = sys.argv[3]
pexp = enfexp[:,list(enfgenes).index(gene)]

# sort columns of individuals in both matrices
indvmask = np.argsort(indv)[np.isin(np.sort(indv),enfindv)]
enfmask = np.argsort(enfindv)[np.isin(np.sort(enfindv),indv)]
pexp, enfindv = pexp[enfmask], enfindv[enfmask]
snp, indv = snp[:, indvmask], indv[indvmask]

# Find drivers that ruin correlation with full model prediction
addname = ''
if '--antidriver' in sys.argv:
    pexp = -pexp
    addname = 'anti'


# check sorting and sizes
print(np.shape(snp), len(avals), len(indv), len(pexp), len(enfindv), np.array_equal(enfindv, indv))

# sort the attributions by size
s = np.argsort(-np.absolute(avals))
avals = np.nan_to_num(avals[s])
avars = avars[s]
snp = snp[s]

# compute percentage of indivuduals with SNV in this set
snpcount = np.sum(snp > 0, axis = 1).astype(int)
snpper = snpcount/len(snp[0])
snpper = np.around(100*snpper,1)

# compute the sum linear approximated prediction from attributions times genotype
imp = snp.T * avals
impm = np.sum(imp, axis = 1)
if '--logsum' in sys.argv:
    impm = np.log10(impm + 1 - np.amin(impm))
# compute the pearson correlation between the linear approximation and the full prediction
peartot = pearsonr(impm, pexp)[0]

print(gene, round(peartot,2))

# Forward selection adds SNVs to sum in order of their absolute attribution
# A driver is an SNV that significantly correlates with the predicted values and also changes the correlation of the partial approximation significantly
if '--forward' in sys.argv:
    obj = open(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversfw.txt','w')
    old_pearsub = 0
    previous = np.zeros(len(impm))
    for i in range(len(s)):
        impn = np.sum(imp[:,:i+1], axis = 1)
        if '--logsum' in sys.argv:
            impn = np.log10(impn + 1-np.amin(impn))
        pearsub = pearsonr(impn, pexp)[0]
        pearexp = pearsub/peartot
        impi = imp[:,i]
        if '--logsum' in sys.argv:
            impi = np.log10(1+impi-np.amin(impi))
        peari = pearsonr(impi, pexp)
        peari, pearisig = peari
        isdriver = False
        if (pearexp - old_pearsub) > 0.05 and peari > 0 and pearisig < 0.01/len(s):
            impj = np.sum(imp[:,np.arange(np.shape(imp)[1]) != i],axis = 1)
            if '--logsum' in sys.argv:
                impj = np.log10(1+impj-np.amin(impj))
            pearj = (peartot - pearsonr(impj, pexp)[0])/peartot
            print(i, str(int(avars[i]))+' '+str(round(avals[i],5))+' '+str(round(peari,2))+' '+ str(round(pearsub,2))+' '+ str(round(pearexp-old_pearsub,2))+' '+ str(snpcount[i])+' '+str(snpper[i]))
            obj.write(str(int(avars[i]))+' '+str(round(avals[i],5))+' '+str(round(peartot,2))+' '+str(round(peari,2))+' '+ str(round(pearsub,2))+' '+str(round(pearj,2))+' '+ str(round(pearexp-old_pearsub,2))+' '+ str(snpcount[i])+' '+str(snpper[i])+'\n')
            isdriver = True
            
        # For understanding: Each selection step can be plotted
        if '--plot_test' in sys.argv:
            tid = int(sys.argv[sys.argv.index('--plot_test')+1])
            if i < tid:
                fig = plt.figure(figsize = (4,4))
                ax = fig.add_subplot(111)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_ylabel('Enformer prediction')
                ax.set_xlabel('Sum prediction')
                ax.set_title(gene+'-'+str(int(avars[i]))+'('+str(i)+')\n'+str(round(avals[i],2))+' Drv=:'+str(isdriver)[0])
                ax.scatter(impm,pexp,color = 'grey', label = 'Full sum '+str(round(peartot,2)))
                ax.scatter(impi,pexp, color = 'goldenrod', label = 'Only SNP '+str(round(peari,2)), alpha = 0.4)
                #ax.scatter(previous, pexp,color='purple', label = 'Sum until SNP-1 '+str(round(old_pearsub,2)), alpha = 0.5)
                ax.scatter(impn, pexp,color='navy', label = 'Sum until SNP '+str(round(pearsub,2)), alpha = 0.5)
                ax.legend()
                fig.savefig(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversfw'+str(i)+'.jpg', dpi = 200, bbox_inches = 'tight')
                print(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversfw'+str(i)+'.jpg')
                plt.close()
        old_pearsub = pearexp
        previous = impn

# Reverse selection removes SNVs form the sum and computes the change of correlation to predicted values
# In every round the SNV with largest reduction of the correlation is selected and the SNV permanently removed from the sum
# In the next round the impact of SNVs is computed on the residual sum without drivers from the previous rounds.
# The process is repeated until the correlation to the predicted values is none-significant
else:
    obj = open(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversbw.txt','w')
    drivers = []
    mask = ~np.isin(np.arange(len(s), dtype = int),drivers)
    pearcomp = pearsonr(impm, pexp)[0]
    pearc = pearsonr(impm, pexp)
    indvimpact = np.ones(len(s))
    indvpear = np.ones(len(s))
    indvpval = np.zeros(len(s))
    removed = 0
    rounds = 0
    while pearc[0] > 0 and pearc[1] < 0.01/len(s):
        pdriver = []
        for i in np.where(mask)[0]:
            if indvpear[i] > 0 and indvpval[i] < 0.01/len(s):
                nmask = np.copy(mask)
                nmask[i] = False
                if rounds == 0:
                    impi = imp[:,i]
                    if '--logsum' in sys.argv:
                        impi = np.log10(1+impi-np.amin(impi))
                    peari = pearsonr(impi, pexp)
                    peari, pearisig = peari
                    indvpear[i] = peari
                    indvpval[i] = pearisig
                    impn = np.sum(imp[:,nmask], axis = 1)
                    if '--logsum' in sys.argv:
                        impn = np.log10(impn + 1-np.amin(impn))
                    pearsub = pearsonr(impn, pexp)[0]
                    indvimpact[i] = 1.-(pearsub/pearcomp)
                    if pearisig < 0.01/len(s) and peari > 0:
                        pdriver.append([i,1.-(pearsub/pearcomp),peari])
                else:
                    impn = np.sum(imp[:,nmask], axis = 1)
                    if '--logsum' in sys.argv:
                        impn = np.log10(impn + 1-np.amin(impn))
                    pearsub = pearsonr(impn, pexp)[0]
                    pdriver.append([i,1.-(pearsub/pearcomp)-removed,indvpear[i]])

        if len(pdriver) == 0:
            break
        pdriver = np.array(pdriver)
        maindriver = np.argmax(pdriver[:,1])
        drivers.append(pdriver[maindriver])
        removed += pdriver[maindriver][1]
        rounds += 1
        mask[int(drivers[-1][0])] = False
        impn = np.sum(imp[:,mask], axis = 1)
        if '--logsum' in sys.argv:
            impn = np.log10(impn + 1-np.amin(impn))
        pearc = pearsonr(impn, pexp)
        print('Round', rounds, pearc, pearcomp)
    
    mask = np.ones(len(s)) == 1
    for dr in drivers:
        i, pearsub, peari = dr
        i = int(i)
        print(i,str(int(avars[i]))+' '+str(round(avals[i],5))+' '+str(round(peari,2))+' '+ str(round(pearsub,2))+' '+ str(round(indvimpact[i],2))+' '+ str(snpcount[i])+' '+str(snpper[i])) 
        obj.write(str(int(avars[i]))+' '+str(round(avals[i],5))+' '+str(round(peari,2))+' '+ str(round(pearsub,2))+' '+ str(round(indvimpact[i],2))+' '+ str(snpcount[i])+' '+str(snpper[i])+'\n')
        # This process can visualized for every round
        if '--plot_test' in sys.argv:
            print(i)
            mask[i] = False
            impn = np.sum(imp[:,mask], axis = 1)
            if '--logsum' in sys.argv:
                impn = np.log10(impn + 1-np.amin(impn))
            fig = plt.figure(figsize = (4,4))
            ax = fig.add_subplot(111)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylabel('Enformer prediction')
            ax.set_xlabel('Sum prediction')
            ax.set_title(gene+'-'+str(int(avars[i]))+'('+str(i)+')\n'+str(round(avals[i],2))+' Drv='+str(isdriver)[0])
            ax.scatter(impm, pexp,color = 'grey', label = 'Full sum '+str(round(peartot,2)))
            ax.scatter(impn, pexp,color='navy', label = 'Sum without SNP '+str(round(pearsub,2)))
            ax.scatter(impi,pexp, color = 'goldenrod', label = 'Only SNP '+str(round(peari,2)))
            ax.legend()
            fig.savefig(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversbw'+str(i)+'.jpg', dpi = 150)
            print(os.path.splitext(sys.argv[1])[0]+'_'+addname+'driversbw'+str(i)+'.jpg')
            plt.close()


