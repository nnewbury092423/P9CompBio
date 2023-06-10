#!/usr/bin/env python3

import numpy as np
import scipy as sp

ERROR_GTR_PARAMS_FILE = "Invalid GTR parameters file"
ERROR_MSA = "Invalid multiple sequence alignment"
ERROR_SCORE = "Invalid log-likelihood score. Must be a negative float"

def likelihood(tree, seqs, gtr_probs, gtr_rates):
    '''
    This function uses Felsenstein's tree-pruning algorithm to compute the log-likelihood score of a tree given an MSA
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param seqs: A dictionary where keys are labels corresponding to the labels of the leaves in ``tree`` and values are sequences (strings)
    :param gtr_probs: The GTR stationary probabilities as a list [prob_A, prob_C, prob_G, prob_T]
    :param gtr_rates: The GTR transition rates as a list [rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]
    :return: The log-likelihood score
    '''
    prob_A, prob_C, prob_G, prob_T = gtr_probs # You can use these if it's more convenient
    rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG = gtr_rates # You can use these if it's more convenient
    log_likelihood_score = 0. # You can modify this variable, or you can use your own; either way is fine
    # TODO Your code here
    # generate Rate Matrix


    AA = -(prob_G*rate_AG+ prob_C*rate_AC + prob_T*rate_AT)
    CC  = -(prob_A*rate_AC+ prob_G*rate_CG + prob_T*rate_CT)
    GG  = -(prob_A*rate_AG+ prob_C*rate_CG + prob_T*rate_GT)
    TT = -(prob_A*rate_AT+ prob_G*rate_GT + prob_C*rate_CT)
    R = np.array([[AA,prob_C*rate_AC, prob_G*rate_AG, prob_T*rate_AT],\
                  [prob_A*rate_AC,CC, prob_G*rate_CG, prob_T*rate_CT],\
                    [prob_A*rate_AG, prob_C*rate_CG,GG, prob_T*rate_GT],\
                    [prob_A*rate_AT, prob_C*rate_CT, prob_G*rate_GT,TT]])

    
    #finaldict
    # normalize R
    #ACGT

    d = 0
    for i in range(len(gtr_probs)):
         d = d - gtr_probs[i] * R[i][i]
    row = 0
    Rnorm = R/abs(d)
    #import pdb; pdb.set_trace()


    L = dict()
    # Recursive Formula
    for node in tree.traverse_postorder():
        # Iniitialize the leaves
        if node.is_leaf():
            store = np.zeros((4,len(seqs[node.get_label()])))
            # Leaf sequences and probabilities are known, so we just generate matrix
            for i,letter in enumerate(seqs[node.get_label()]):
                if letter == 'A':
                    store[:,i] = [1,0,0,0]
                elif letter == 'C':
                    store[:,i] = [0,1,0,0]
                elif letter =='G':
                    store[:,i] = [0,0,1,0]
                elif letter == 'T':
                    store[:,i] = [0,0,0,1]
            L[node] = store


    #Main recursion
        else:
            prod = 1
            #Multiply each child node likelihood matrix together elementwise because each event is independent
            for child in node.child_nodes():
                #calculate probability matrices
                t = child.get_edge_length() 
                tranP = sp.linalg.expm(t*Rnorm)
                # find the matrix from the child node
                childprobs = L[child]

                # tranP * childprobs matrix multiplication fills out the elements
                prod = prod *np.matmul(tranP,childprobs)
            L[node] = prod
    
    # for the root, we do the same with the independent probabilities pi
    Likear = np.matmul(gtr_probs,L[tree.root])

    # to find the total log likelihood, we take the log sum
    log_likelihood_score = np.sum(np.log(Likear))
            
    return log_likelihood_score


def returnindex(letter):
    if letter == 'A':
        index = 0
    elif letter == 'C':
        index = 1
    elif letter =='G':
        index = 2
    elif letter == 'T':
        index = 3


    return index

def read_FASTA(filename):
    '''
    This function reads a FASTA file from a file and returns a dictionary mapping identifiers to sequences
    '''
    stream = open(filename); seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq; stream.close()
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your functions, and outputting the results
    '''
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-p', '--gtr_params', required=True, type=str, help="GTR Parameters File")
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Multiple Sequence Alignment (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Log-Likelihood Score")
    args = parser.parse_args()

    # load input tree
    from treeswift import read_tree_newick
    if args.tree == 'stdin':
        from sys import stdin
        tree = read_tree_newick(stdin)
    else:
        tree = read_tree_newick(args.tree)

    # load GTR parameters
    try:
        gtr_probs, gtr_rates = [[float(e.strip()) for e in l.strip().split()] for l in open(args.gtr_params)]
    except:
        raise ValueError(ERROR_GTR_PARAMS_FILE)

    # load MSA
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_MSA)
    for l in tree.traverse_leaves():
        assert l.label in seqs, "Missing sequence for: %s" % l.label

    # run student code and output
    score = likelihood(tree, seqs, gtr_probs, gtr_rates)
    if (not isinstance(score,float) and not isinstance(score,int)) or score > 0:
        raise ValueError(ERROR_SCORE)
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    outfile.write('%s\n' % str(score)); outfile.close()
