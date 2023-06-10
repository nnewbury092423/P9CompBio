
import numpy as np
import compute_likelihood as cp
import os
import joblib
from treeswift import read_tree_newick
import random
from operator import itemgetter
import csv
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib.colors as mcolors
mainpath = 'one'
num = 5000


def resample(seqs):
    seqsnew = dict()
    keyslist = list(seqs.keys())
    
    filterindices = np.random.choice(len(seqs[keyslist[0]]),len(seqs[keyslist[0]]))
    for key in seqs:
        
        alnstring = seqs[key]  
        seqsnew[key] = ''.join(itemgetter(*filterindices)(alnstring))

    return seqsnew

def plothistogram(likelihoods,Ln):
    #print(likelihoods)
    colorlist= list(mcolors.CSS4_COLORS.keys())  
    i = 0
    for k1 in likelihoods:
        #print(likelihoods[k1])
        plt.axvline(x=likelihoods[k1], label = k1, color = colorlist[i])
        #plt.text(10.1,0,k1,rotation=90)
        i+=2

    counts, bins = np.histogram(Ln)
    plt.stairs(counts, bins)
    plt.xlim([-36000,-33000])
    #plt.show(block=True)
    plt.legend(loc = 'upper right',fontsize = "5")
    plt.savefig("test.png")




if __name__ == "__main__":
#generate a bunch of sequences
#throw a bunch at likelihood
#get main param
    tot = [None]*10

    for senum in range(1,11):
        if senum <10:
            seqs = cp.read_FASTA('example/0'+str(senum) + '.fas')
        if senum == 10:
            seqs = cp.read_FASTA('example/10.fas')

    #get main tree
        maintree = read_tree_newick('example/tree.nwk')
        gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open('example/gtr_params.txt')]
        #import pdb; pdb.set_trace()
        
        # compute the null distribution
        Ln = [None] *num
        for i in range(num):
            seqsnew = resample(seqs)
            #import pdb; pdb.set_trace()
            Ln[i] = cp.likelihood(maintree,seqsnew, gtr_probs, gtr_rates)

        #import pdb; pdb.set_trace()

        #sort

        Ln.sort()
        La = dict()
        likelihoods = dict()
        # for main
        La['gtr_params.tx'] = {}
        La['gtr_params.tx']['t tree.nw'] = np.searchsorted(Ln,cp.likelihood(maintree, seqs, gtr_probs, gtr_rates))/num
        likelihoods['gtr_params.txt tree.nw'] = cp.likelihood(maintree, seqs, gtr_probs, gtr_rates)
        #import pdb; pdb.set_trace()  
    #loop through params with main tree
        
        for filename in os.listdir('params'):
            f = os.path.join('params', filename)
            # checking if it is a file
            if os.path.isfile(f):
                gtr_probsnew, gtr_ratesnew = [[float(e.strip()) for e in l.strip().split()] for l in open(f)]
                La[filename] = {} 
                # a 1 means that we found better likelihoods
                like= cp.likelihood(maintree, seqs, gtr_probsnew, gtr_ratesnew)
                La[filename][' tree.nwk']=  np.searchsorted(Ln,like)/num
                likelihoods[filename + 'tree.nwk'] = like       
        La['gtr_params.txt '] = {}
        for filename in os.listdir('trees'):
            f = os.path.join('trees', filename)
            # checking if it is a file
            if os.path.isfile(f):
                tree = read_tree_newick(f)
                like2 = cp.likelihood(tree, seqs, gtr_probs, gtr_rates)
                La['gtr_params.txt '][filename] = np.searchsorted(Ln,like2)/num
                likelihoods['gtr_params.txt' + filename] =  like2


        #import pdb; pdb.set_trace()

        tot[senum-1] = La 
        print(senum)
    # loop through main tree with params    


    import pdb; pdb.set_trace()
    plothistogram(likelihoods,Ln)
    
    with open('output-example.csv', 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for i in range(len(tot)):
            elem = tot[i]
            for key1 in elem:
                for key2 in elem[key1]:
                    stringtowrite = str(i+1)+' '+ key1 + key2+ ' '+str(elem[key1][key2])
                    spamwriter.writerow([stringtowrite])
    
                
   
        
    
    import pdb; pdb.set_trace()





# save to csv

#loop through trees with main params





    #print('hello')