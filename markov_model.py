import pdb
import math
import re

def calc_mm_score(x, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs):
    log_p_x = calc_prob(x, k, trusted_genes, trusted_genes_probs)
    log_q_x = calc_prob(x, k, trusted_non_genes, trusted_non_genes_probs)
    return log_p_x - log_q_x

# k = 3 --> P(x) = P(x1 x2 x3) * P(x4 | x1 x2 x3) * P(x5 | x2 x3 x4) ... P(xn | xn-3 xn-2 xn-1)
def calc_prob(x, k, training, probs):
    prob = probs[x[:k]]
    log_prob = math.log(prob, 2)
    for idx in xrange(0, len(x)-k):
        sub_prob = probs[x[idx:idx+k+1]]
        log_prob += math.log(sub_prob, 2)
    return log_prob

# k = 3 --> P(x) = P(AAA), e.g.
# assume x is a sequence of k nucleotides
def uncondit_prob(x, k, training):
    num_x = 0
    total = 0
    for seq in training:
        for i in xrange(len(seq)-(k-1)):
            sub_seq = seq[i:i+k]
            if sub_seq == x:
                num_x += 1
            total += 1
    return num_x/(total * 1.0)

def condit_prob(x, y, k, training):
    # P(x|y) = P(xy)/P(y)
    y_seqs = find_y(y, k, training)
    num_y = len(y_seqs)
    num_yx = count_yx(y+x, k, y_seqs)
    # P(xy) = Count(xy)/Count(total), P(y) = Count(y)/Count(total) --> P(xy)/P(y) = Count(xy)/Count(y)
    return num_yx/(num_y*1.0)

# find all seqs y1y2y3_
def find_y(y, k, training):
    y_seqs = []  
    [find_y_helper(y, k, seq, y_seqs) for seq in training]  
    return y_seqs    

def find_y_helper(y, k, seq, y_seqs):
    for i in xrange(len(seq)-k):
        sub_seq = seq[i:i+k+1]
        if sub_seq[:k]==y:
            y_seqs.append(sub_seq)

# count seqs y1y2y3x1
def count_yx(yx, k, seqs): 
    count = 0 
    for seq in seqs: 
        for i in xrange(len(seq)-k):
            sub_seq = seq[i:i+k+1]
            if sub_seq==yx: 
                count += 1
    return count     

def calc_probs(k, training):
    NUCLEOTIDES = ["A", "C", "G", "T"]
    probs = {}
    for nuc1 in NUCLEOTIDES: 
        for nuc2 in NUCLEOTIDES:
            for nuc3 in NUCLEOTIDES: 
                y = nuc1 + nuc2 + nuc3 
                probs[y] = uncondit_prob(y, k, training)
                for nuc4 in NUCLEOTIDES:
                    x = nuc4
                    yx = y+x 
                    probs[yx] = condit_prob(x, y, k, training)
    return probs