# Align sequences and cluster them to find common motifs in sequences
# python3 ../find_common_motifs.py ALLgenes_ism_attributions_driversfw_refseq_winsize13.npz complete 0.25 5 1 --geneset ../Prediction_correlationsCageAdultBrainGeneSpecific_CorrelationtoObsRandomNull_tstats_set0.2.list

import numpy as np
import sys, os
from sklearn.cluster import AgglomerativeClustering
    
# function that iterates over all different shifts of two sequences to each other and determines which stretches fullfil minsize, maxgap requirements, scores them and then returns list with valid alignments
def find_stretches(mseq, minsize, maxgap):
    stretches = []
    stretchscores = []
    # stretches need a minimal size and maximum gap and at least a coverge resulting from both
    mincover = (minsize-maxgap)/minsize
    # iterate over all possible start positionts
    mstr = 0
    for l in range(len(mseq)-minsize):
        score = 0
        gap = 0
        #check = np.zeros(len(mseq)-l)
        if mseq[l] == 1:
            for k in range(l, len(mseq)):
                if mseq[k] > 0:
                    gap = 0
                    score += mseq[k]
                    #check[k-l] = 1
                else:
                    gap += 1
                if gap > maxgap:
                    # stop iteration of maxgap is reached
                    break
            # chck if minimum length and coverage is given
            lenstr = k - gap -l +1
            if lenstr >= minsize and score/(lenstr) >= mincover and ((k-gap) > mstr):
                #print(score, gap, lenstr, check)
                stretches.append(np.arange(l,k - gap+1, dtype = int))
                stretchscores.append(score)
                mstr = k-gap
    return stretches, np.array(stretchscores)

# Compute the best alignment score and the location of the matching bases in both sequences between all sequences in seqs
# minsize : minimum size of matching stretches with maxgap positions 
# maxgap : maximum tolerated gaps in stretch
# always_include : determines if the center should alwasy be included in the aligned stretch
def alignNscore_seqpairs(seqs, minsize = 5, maxgap =1, always_include = 'center'):
    nseqs = len(seqs)
    lseqs = len(seqs[0])
    
    # decide if specific position always has to be in the aligned regions
    if always_include == 'center':
        always_include = int(lseqs/2)
    elif always_include == 'none':
        always_include = False
    
    # Generate a matrix with the alignment scores 
    corrmat = np.zeros((nseqs, nseqs))
    # save the stretches that these scores represent
    # dimensions of strmat: (two_compared_seqs, nseqs, nseqs, start_and_end)
    strmat = np.zeros((2,nseqs,nseqs,2), dtype = int)
    for s, seq in enumerate(seqs):
        if s % 50 == 0:
            print(s)
        for q, qeq in enumerate(seqs[s+1:]):
            t = q + s +1
            # iterate over different shifts of the sequences to each other
            stretch0 = [0]
            stretch1 = [0]
            strscore = 0
            for l in range(lseqs - minsize+1):
                # compute positional overlap
                mseq = np.sum(seq[l:] * qeq[:lseqs-l], axis = 1)
                # determine longest aligned stretch with maximal gap of maxgap
                potstretches, potstretchscore = find_stretches(mseq, minsize, maxgap)
                # if a position should always be included check if they are
                if len(potstretchscore) > 0:
                    if always_include != 'none':
                        posinc = np.unique([always_include-l, always_include])
                        keep = [np.sum(np.isin(pot, posinc))==len(posinc) for pot in potstretches]
                        potstretches = [potstretches[p] for p in np.where(keep)[0]]
                        potstretchscore = potstretchscore[keep]
                    # choose longest stretch with mostoverlap and save it
                    if len(potstretchscore) > 0:
                        bscore = np.argmax(potstretchscore)
                        if potstretchscore[bscore] > strscore:
                            # maintain best score, and stretches for both sequences
                            strscore = potstretchscore[bscore]
                            stretch0 = potstretches[bscore]+l
                            stretch1 = potstretches[bscore]
                
                # and repeat with shifts in other direction
                mseq = np.sum(seq[:lseqs-l] * qeq[l:], axis = 1)
                # determine longest aligned stretch with maximal gap of maxgap
                potstretches, potstretchscore = find_stretches(mseq, minsize, maxgap)
                # if a position should always be included check if they are
                if len(potstretchscore) > 0:
                    if always_include != 'none':
                        posinc = np.unique([always_include-l, always_include])
                        keep = [np.sum(np.isin(pot, posinc)) == len(posinc) for pot in potstretches]
                        potstretches = [potstretches[p] for p in np.where(keep)[0]]
                        potstretchscore = potstretchscore[keep]
                        
                    # choose longest stretch with mostoverlap and save it
                    if len(potstretchscore) > 0:
                        bscore = np.argmax(potstretchscore)
                        if potstretchscore[bscore] > strscore:
                            strscore = potstretchscore[bscore]
                            stretch0 = potstretches[bscore]
                            stretch1 = potstretches[bscore]+l
            corrmat[s,t] = corrmat[t,s] = strscore
            strmat[0,s,t] = strmat[1,t,s] = [np.amin(stretch0), np.amax(stretch0)+1]
            strmat[1,s,t] = strmat[0,t,s] = [np.amin(stretch1), np.amax(stretch1)+1]
            #print(s,t,corrmat[s,t], stretch0, stretch1, strmat[0,s,t], strmat[1,s,t])
    return corrmat, strmat

# generate pwms from the assigned clusters and the seqs
# seqs : one hot encoded sequences
# clusters : cluster that seqs belong to
# scoremat : cluster score between sequences to select the most similar sequences as seed
# stretches : positions in the sequence that are aligned and fulfill the requirements
# find : determine where in the pwm the driver snp is located
def generate_pwms(seqs, clusters, scoremat, stretches, find = None):
    clusterpwms = []
    clusteralignments = []
    nts = np.array(list('ACGT'))
    lseqs, nb = np.shape(seqs[0])
    lseqoff = int((lseqs-1)/2)
    print(lseqoff)
    if find is not None:
        findloc = -np.ones(len(seqs), dtype = int)
    for c in np.unique(clusters):
        mask = np.where(clusters == c)[0]
        #print(c, mask)
        if len(mask) > 1:
            seed = mask[np.argmax(np.sum(scoremat[mask][:,mask],axis =1))]
            #print(seed)
            maskorder = mask[np.argsort(-scoremat[seed,mask])]
            #print(maskorder, -scoremat[seed,maskorder])
            fst, fen = stretches[0,seed,maskorder[0],0], stretches[0,seed,maskorder[0],1]
            if find is not None:
                findloc[seed] = lseqoff + find
            pwm = np.zeros((lseqs+2*lseqoff,nb))
            align = np.chararray((len(mask), lseqs+2*lseqoff), itemsize = 1, unicode = True)
            align[:] = '-'
            pwm[fst+lseqoff:fen+lseqoff] += seqs[seed][fst:fen]
            align[0,lseqoff:lseqoff+lseqs] = nts[np.where(seqs[seed]==1)[1]]
            ji = 0
            for m in maskorder:
                if m != seed:
                    ji +=1
                    fst, fen = stretches[0,seed,m,0], stretches[0,seed,m,1]
                    st, en = stretches[1,seed,m,0], stretches[1,seed,m,1]
                    #print(seqs[seed][fst:fen])
                    #print(seqs[m][st:en])
                    #print(pwm[fst+lseqoff:fen+lseqoff])
                    pwm[fst+lseqoff:fen+lseqoff] += seqs[m][st:en]
                    align[ji,lseqoff+fst-st:lseqoff+lseqs+fst-st] = nts[np.where(seqs[m]==1)[1]]
                    #print(pwm[fst+lseqoff:fen+lseqoff])
                    if find is not None:
                        findloc[m] = lseqoff + find + fst - st
            pmask = np.sum(pwm,axis=1)>0
            if find is not None:
                findloc[mask]-=np.where(pmask)[0][0]
            pwm = pwm[np.sum(pwm,axis=1)>0]
            pwm = pwm/np.amax(np.sum(pwm,axis = 1))
        else:
            pwm = None
            align = [nts[np.where(seqs[mask[0]]==1)[1]]]
        clusterpwms.append(pwm)
        clusteralignments.append(align)
    return clusterpwms, clusteralignments, findloc

# write pwm wiht pwmname and alphabet into output_file_path
def write_meme_file(pwm, pwmname, alphabet, output_file_path):
    """[summary]
    write the pwm to a meme file
    Args:
        pwm ([np.array]): n_filters * 4 * motif_length
        output_file_path ([type]): [description]
    """
    n_filters = len(pwm)
    meme_file = open(output_file_path, "w")
    meme_file.write("MEME version 4 \n")
    meme_file.write("ALPHABET= "+alphabet+" \n")
    meme_file.write("strands: + -\n")

    print("Saved PWM File as : {}".format(output_file_path))

    for i in range(0, n_filters):
        if pwm[i] is not None:
            if np.sum(pwm[i]) > 0:
                meme_file.write("\n")
                meme_file.write("MOTIF %s \n" % pwmname[i])
                meme_file.write("letter-probability matrix: alength= "+str(len(alphabet))+" w= %d \n"
                        % np.count_nonzero(np.sum(pwm[i], axis=0)))
        
            for j in range(0, np.shape(pwm[i].T)[-1]):
                for a in range(len(alphabet)):
                    if a < len(alphabet)-1:
                        meme_file.write(str(pwm[i].T[ a, j])+ "\t")
                    else:
                        meme_file.write(str(pwm[i].T[ a, j])+ "\n")

    meme_file.close()


if __name__ == '__main__':
    
    # read one-hot encoded seqfile as npz
    seqfile = sys.argv[1]
    
    seqf = np.load(seqfile)
    seqs = seqf['seqs']
    genes = seqf['genes']
    loci = seqf['loci']
    print(np.shape(seqs))

    # define linkage for agglomerative clustering
    linkage = sys.argv[2]
    # define cutoff for clusters
    cutoff = float(sys.argv[3])
    # define minimum size of aligned stretch of sequences
    msize = int(sys.argv[4])
    # define maximum allowed gap in stretch of aligned sequences
    maxgap = int(sys.argv[5])

    # select subset of genes
    if '--geneset' in sys.argv:
        geneset = np.genfromtxt(sys.argv[sys.argv.index('--geneset')+1], dtype = str)
        gmask = np.isin(genes, geneset)
        seqs = seqs[gmask]
        genes = genes[gmask]
        loci = loci[gmask]

    outname = os.path.splitext(seqfile)[0]+'_clust_ms'+str(msize)+'-'+str(maxgap)+'_'+linkage+str(cutoff)

    # generate alignment and return scoring matrix between all sequences and a matrix with the location of the stretches that are aligned
    scoremat, alignmat = alignNscore_seqpairs(seqs, minsize = msize, maxgap =maxgap, always_include = 'center')
    # generate distance matrix from score matrix
    # combine everything that has more or as many as msize similarity
    distmat = 1. - scoremat/msize
    distmat[distmat < 0] = 0
    np.fill_diagonal(distmat, 0)

    # apply agglomerative clustering
    clustering = AgglomerativeClustering(n_clusters = None, affinity = 'precomputed', linkage = linkage, distance_threshold = cutoff).fit(distmat)
    clusters = clustering.labels_
    
    # save assigned clusterf os sequences
    np.savetxt(outname + '_clusters.txt', np.array([genes, loci, clusters]).T, fmt = '%s')
    print(len(seqs), 'form', len(np.unique(clusters)), 'clusters')
    
    # generate pwms for clusters and determine the location of the center in these pwms
    clusterpwms, clusteralignments, centerloc = generate_pwms(seqs, clusters, scoremat, alignmat, find = int(len(seqs[0])/2))
    # save location of the snp in the pwms
    np.savetxt(outname + '_locclpwms.txt', np.array([genes, loci, centerloc]).T, fmt = '%s')
    clusternames = ["Cluster_"+str(i) for i in range(len(clusterpwms))]
    # save pwms in meme format
    write_meme_file(clusterpwms, clusternames, 'ACGT', outname+'_clusterpwms.txt')
    # print stats
    for c in np.unique(clusters):
        print('Cluster', c, int(np.sum(clusters == c)))
        for cl in clusteralignments[c]:
            print( ''.join(cl))
        #print(scoremat[clusters == c][:, clusters == c])



