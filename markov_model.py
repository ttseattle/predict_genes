import pdb
import math

NUCLEOTIDES = ["A", "C", "G", "T"]

def calc_mm_score(x, k, trusted_genes, trusted_non_genes):
    p_x, log_p_x = calc_prob(x, k, trusted_genes)
    q_x, log_q_x = calc_prob(x, k, trusted_non_genes)
    mm_score = log_p_x - log_q_x    
    return mm_score

# TODO: Exclude the stop codon from your frequency estimates in both "trusted" sets; it might introduce bias, especially in the short negative examples.
# k = 3 --> P(x) = P(x1 x2 x3) * P(x4 | x1 x2 x3) * P(x5 | x2 x3 x4) ... P(xn | xn-3 xn-2 xn-1)
def calc_prob(x, k, training):
    prob = uncondit_prob(x[:k], k, training)
    log_prob = math.log(prob, 2)

    idx = 0 
    while idx < len(x)-k:
        sub_prob = condit_prob(x[idx+k], x[idx:idx+k], k, training)
        prob *= sub_prob 
        log_prob += math.log(sub_prob, 2)
        idx += 1

    return (prob, log_prob)

# k = 3 --> P(x) = P(AAA), e.g.
# assume x is a sequence of k nucleotides 
def uncondit_prob(x, k, training):
    num_x = 0 
    total = 0 
    for orf in training:
        seq = orf[0]
        for i in range(len(seq)-(k-1)):
            sub_seq = seq[i:i+k]
            if sub_seq == x:
                num_x += 1 
            total += 1 
    return num_x/(total * 1.0)

def condit_prob(x, y, k, training):
    y_seqs = find_y(y, k, training)
    num_y = len(y_seqs)
    num_x = count_x(x, k, y_seqs)
    return num_x/(num_y * 1.0)

# find all seqs y1y2y3_ 
def find_y(y, k, training):
    y_seqs = []
    for orf in training:
        seq = orf[0]
        for i in range(len(seq)-k):
            sub_seq = seq[i:i+k+1]
            if sub_seq[:k]==y:
                y_seqs.append(sub_seq)
    return y_seqs 

# find all seqs ___x
def count_x(x, k, seqs):
    num_x = 0 
    for seq in seqs:
        if seq[-1]==x:
            num_x += 1
    return num_x