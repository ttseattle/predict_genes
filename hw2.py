import open_reading_frame 
import markov_model
import matplotlib.pyplot as plt
import numpy as np
import pdb
from numpy import trapz
import time
start_time = time.time()

MJANNASCHII_FILEPATH = "data/genome.fna"
ANNOTATIONS_FILEPATH = "data/annotation.gbff"
SHORT_THRESH = 50
LONG_THRESH = 1400
K = 3 



def main(): 
    # get data 
    seq = open_reading_frame.prep_mjannaschii(MJANNASCHII_FILEPATH) 
    gene_ends = open_reading_frame.prep_annotations(ANNOTATIONS_FILEPATH, seq)

    orfs = find_orfs(seq)
    short_orfs, long_orfs = find_short_and_long_orfs(orfs, SHORT_THRESH, LONG_THRESH)
    get_num_cdss(gene_ends)

    trusted_non_genes = prep_training_set(short_orfs)
    trusted_genes = prep_training_set(long_orfs)
    trusted_genes_probs = markov_model.calc_probs(K, trusted_genes)
    trusted_non_genes_probs = markov_model.calc_probs(K, trusted_non_genes)

    # homework problems 
    

    # scores = orf_mm_scores(orfs, gene_ends, 30, K, trusted_genes, trusted_non_genes)
    # print mm_thresh_2(orfs, gene_ends, 30, K, trusted_genes, trusted_non_genes, scores)
    # print mm_thresh_2(orfs, gene_ends, 40, K, trusted_genes, trusted_non_genes, scores)
    # print mm_thresh_2(orfs, gene_ends, 50, K, trusted_genes, trusted_non_genes, scores)
    

    # open_reading_frame.len_thresh(orfs, gene_ends, 348)
    # 
    # count_xy(K, trusted_genes, trusted_non_genes)
    print_summary(short_orfs, long_orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)

# 1A #####################################################################

def find_orfs(seq): 
    print "NUMBER OF ORFS FOUND"
    orfs = []

    for offset in xrange(3):
        orfs_at_offset = open_reading_frame.scan(seq, offset)
        orfs += orfs_at_offset
        print "number of ORFS at offset " + str(offset) + ":", len(orfs_at_offset)

    print "total number of ORFS:", len(orfs)
    print ""
    return orfs

# 1B & 1C ################################################################

def find_short_and_long_orfs(orfs, short_threshold, long_threshold):
    print "NUMBER OF ORFS OF LENGTH < 50 AND > 1400"
    trusted_non_genes = []
    trusted_genes = []
    for orf in orfs: 
        seq, start, end = orf
        if len(seq) < short_threshold:
            trusted_non_genes.append(orf)
        elif len(seq) > long_threshold:
            trusted_genes.append(orf)

    print "less than 50:", len(trusted_non_genes)
    print "greater than 1400:", len(trusted_genes)
    print ""
    return (trusted_non_genes, trusted_genes)

# 1D #####################################################################

def get_num_cdss(cdss):
    print "NUMBER OF SIMPLE FORWARD STRAND CDSs"
    print "number of genes in GenBank:", len(cdss)
    print ""

# 1E #####################################################################

def remove_stop_codons(training):
    new_training = []
    for orf in training: 
        seq = orf[0]
        new_seq = seq[:-3]
        new_training.append((new_seq, orf[1], orf[2]))
    return new_training

def remove_stop_codon(orf):
    seq = orf[0]
    return seq[:-3]

def prep_training_set(training): 
    new_training = []
    for orf in training: 
        seq = orf[0]
        new_seq = seq[:-3]
        if len(new_seq) != 0:
            new_training.append(new_seq)
    return new_training

def count_xy(k, trusted_genes, trusted_non_genes): 
    print "COUNTS FOR CALCULATING P(A|Axy) -- COUNT(AxyA), COUNT(Axy)"
    count_xy_helper(k, trusted_genes)

    print "COUNTS FOR CALCULATING Q(A|Axy) -- COUNT(AxyA), COUNT(Axy)"
    count_xy_helper(k, trusted_non_genes)

def count_xy_helper(k, training):
    print "{:2} {:<15} {:<15} {:<15} {:<15}".format("", "A", "C", "G", "T")
    for nuc1 in ["A", "C", "G", "T"]:
        counts = []
        for nuc2 in ["A", "C", "G", "T"]:
            axy = 'A' + nuc1 + nuc2 
            y_seqs = markov_model.find_y(axy, k, training)
            num_y = len(y_seqs)
            num_x = markov_model.count_x('A', k, y_seqs)
            counts.append(str(num_x) + ", " + str(num_y))
        print "{:2} {:<15} {:<15} {:<15} {:<15}".format(nuc1, counts[0], counts[1], counts[2], counts[3])
    print ""

# 1F #####################################################################

def print_summary(short_orfs, long_orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends):
    print "SUMMARY FOR FIRST 5 ORFS OF LENGTH < " + str(SHORT_THRESH)
    print_summary_helper(short_orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)
    print ""

    print "SUMMARY FOR FIRST 5 ORFS OF LENGTH > " + str(LONG_THRESH)
    print_summary_helper(long_orfs, trusted_genes, trusted_non_genes,trusted_genes_probs, trusted_non_genes_probs, gene_ends)
    print ""

def print_summary_helper(orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends):
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format("#", "Start", "Length", "MM Score", "CDS?")
    lens = []
    mms = []
    count = 1
    for orf in sorted(orfs, key=lambda tup: tup[1])[:5]:
        len, mm = print_info(count, orf, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)
        lens.append(len)
        mms.append(mm)
        count += 1

def print_info(count, orf, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends): 
    seq, start, end = orf 
    mm_score = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
    cds = 'True' if end in gene_ends else 'False' 
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format(count, start, len(seq), mm_score, cds)
    return (len(seq), mm_score)

# 2 ######################################################################

def plot_roc():
    len_x, len_y = plot_len_roc(1400, orfs, gene_ends)
    plt.plot(len_x, len_y, 'ro')
    
    mm_x, mm_y = plot_mm_roc(100, K, orfs, gene_ends, trusted_genes, trusted_non_genes)
    plt.plot(mm_x, mm_y, 'ro')
    
    # Compute the area using the composite trapezoidal rule.
    len_area = trapz(len_y, dx=5)
    print("len ROC area under curve =", len_area)
    mm_area = trapz(mm_y, dx=5)
    print("mm ROC area under curve =", mm_area)
    
    plt.show()

def orf_mm_scores(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes):
    mm_scores = []
    for orf in orfs: 
        seq, start, end = orf 
        mm_score = markov_model.calc_mm_score(seq, k, trusted_genes, trusted_non_genes)
        mm_scores.append((orf, mm_score))
    return mm_scores 

def plot_len_roc(threshold, orfs, gene_ends):
    x = []
    y = []
    thresh = 0 
    while thresh < threshold:
        thresh += 50
        x_pt, y_pt = len_thresh(orfs, gene_ends, thresh)
        x.append(x_pt)
        y.append(y_pt)
    return x, y 

def plot_mm_roc(threshold, k, orfs, gene_ends, trusted_genes, trusted_non_genes): 
    x = []
    y = []
    thresh = 0 
    while thresh < threshold: 
        thresh += 50
        x_pt, y_pt = mm_thresh(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes)
        x.append(x_pt)
        y.append(y_pt)
    return x, y 

def get_rates(num_true_positives, num_false_positives, orfs, genes, gene_ends): 
    # TPR = TP/TP+FN 
    # TP + FN = CDS 
    tp_rate = num_true_positives*1.0/len(gene_ends)

    # FPR += FP/FP+TN 
    # TN = orfs not meeting threshold or in cds 
    fp_rate = num_false_positives*1.0/(len(orfs)-len(genes)-len(gene_ends)+num_true_positives)
    
    return (fp_rate, tp_rate)

# 3 ######################################################################

# true-positive: a gene, and says is a gene
def positive(orfs, gene_ends):
    num_false_positives = 0
    num_true_positives = 0
    for orf in orfs:
        seq, start, end = orf 
        if end in gene_ends:
            num_true_positives += 1
        else: 
            num_false_positives +=1
    return (num_false_positives, num_true_positives)

def len_thresh(orfs, gene_ends, threshold):
    genes = []
    for orf in orfs:
        seq, start, end = orf 
        if len(seq) > threshold:
            genes.append(orf)
    
    num_false_positives, num_true_positives = positive(genes, gene_ends)

    return get_rates(num_true_positives, num_false_positives, orfs, genes, gene_ends)
    
# 4 ######################################################################
 
 def mm_thresh(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes):
    genes = []
    for orf in orfs: 
        seq, start, end = orf 
        mm_score = markov_model.calc_mm_score(seq, k, trusted_genes, trusted_non_genes)
        if mm_score > threshold: 
            genes.append(orf)

    num_false_positives, num_true_positives = positive(genes, gene_ends)

    return get_rates(num_true_positives, num_false_positives, orfs, genes, gene_ends)

def mm_thresh_2(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes, mm_scores):
    genes = []
    for mm_score in mm_scores: 
        orf, score = mm_score 
        if score > threshold: 
            genes.append(orf) 
    num_false_positives, num_true_positives = positive(genes, gene_ends)
    return get_rates(num_true_positives, num_false_positives, orfs, genes, gene_ends) 

# ERROR CHECKING #########################################################

def check_genes(orfs, gene_ends):
    count = 0
    for orf in orfs: 
        subseq, start, end = orf
        if end in gene_ends:
            count += 1
    print "matching ends:", count 
    print ""

# EXECUTE PROGRAM ########################################################

main()
print("--- %s seconds ---" % (time.time() - start_time))

# x1 = [39, 6, 15, 12, 15]
# x2 = [1614, 1623, 2388, 1422, 1665]
# y1 = [1.56678727191, 0.667366861992, 0.633938280251, 0.418546201372, 0.634136080416]
# y2 = [54.5409415958, 65.6859640015, 76.0581583166, 36.7210402482, 61.6096434786]
# plt.scatter(x1, y1, color='blue')
# plt.scatter(x2, y2, color='orange')

# x1_med = x1[2]
# y1_med = y1[2]
# x2_med = x2[2]
# y2_med = y2[2]
# print "(" + str(x1_med) + ", " + str(y1_med) + ")"
# print "(" + str(x2_med) + ", " + str(y2_med) + ")"

# plt.plot([x1_med, x2_med], [y1_med, y2_med])

# plt.show()

# x = [0.56, 0.78, 0.91]
# y = [0.01, 0.19, 0.58]
# plt.plot(y, x, 'ro')
# plt.show()

