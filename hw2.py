import open_reading_frame
import markov_model
import matplotlib.pyplot as plt
import pdb
import numpy as np
import time

MJANNASCHII_FILEPATH = "data/genome.fna"
ANNOTATIONS_FILEPATH = "data/annotation.gbff"
SHORT_THRESH = 50
LONG_THRESH = 1400
K = 3

def main():
    # get data
    seq = open_reading_frame.prep_mjannaschii(MJANNASCHII_FILEPATH)
    gene_ends = open_reading_frame.prep_annotations(ANNOTATIONS_FILEPATH, seq)

    # get ORFs
    orfs = find_orfs(seq)

    # get training data
    short_orfs, long_orfs = find_short_and_long_orfs(orfs, SHORT_THRESH, LONG_THRESH)
    get_num_cdss(gene_ends)
    trusted_non_genes = remove_stop_codons(short_orfs)
    trusted_genes = remove_stop_codons(long_orfs)

    # get probabilities for calculating Markov Model scores
    trusted_genes_probs = markov_model.calc_probs(K, trusted_genes)
    trusted_non_genes_probs = markov_model.calc_probs(K, trusted_non_genes)

    # get information about 5 short and 5 long ORFs
    # count_xy(K, trusted_genes, trusted_non_genes)
    # print_summary(short_orfs, long_orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)

    # # find thresholds to reach 80% true positive rate
    # reach_len_thresh(orfs, gene_ends, 434) # 435 --> 79.72665148063781%
    # reach_mm_thresh(orfs, gene_ends, 6.774, K, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs) # 6.775 --> 79.95444191343963%

    # # plot data
    plot_roc(orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)

    # plot_scatter(K, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends, 0.20)

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

def plot_roc(orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends):
    plt.suptitle('ROC Curves')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    axes = plt.gca()
    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])

    len_x, len_y = plot_len_roc(orfs, gene_ends)
    plt.plot(len_x, len_y, color='red')

    mm_x, mm_y = plot_mm_roc(K, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)
    plt.plot(mm_x, mm_y, color='green')

    combo_x, combo_y = plot_combo_roc(K, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)
    plt.plot(combo_x, combo_y, color='magenta')

    # Compute the area using the composite trapezoidal rule.
    len_area = np.trapz(len_y, sorted(len_x))
    print "len ROC area under curve =", len_area
    mm_area = np.trapz(sorted(mm_y), sorted(mm_x))
    print "mm ROC area under curve =", mm_area
    combo_area = np.trapz(sorted(combo_y), sorted(combo_x))
    print "combo ROC area under curve =", combo_area

    plt.legend(['length', 'markov model', 'combined length & markov model'], loc='lower right')
    plt.show()

def plot_len_roc(orfs, gene_ends):
    x = []
    y = []
    thresh = 0
    while thresh < LONG_THRESH:
        thresh += 10
        data = len_thresh(orfs, gene_ends, thresh)
        x_pt = data[2]
        y_pt = data[3]
        x.append(x_pt)
        y.append(y_pt)
    return x, y

def plot_mm_roc(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends):
    x = []
    y = []
    thresh = 0
    while thresh < 20:
        thresh += 0.2
        data = mm_thresh(orfs, gene_ends, thresh, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
        x_pt = data[2]
        y_pt = data[3]
        x.append(x_pt)
        y.append(y_pt)
    return x, y

def plot_combo_roc(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends): 
    trusted_gene_x, trusted_gene_y = plot_scatter_helper(k, trusted_genes, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)

    trusted_non_gene_x, trusted_non_gene_y = plot_scatter_helper(k, trusted_non_genes,trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)

    a_x, a_y, b_x, b_y, m, b = calc_median(trusted_gene_x, trusted_gene_y, trusted_non_gene_x, trusted_non_gene_y)
    
    x = []
    y = []
    thresh = 0
    while thresh < 1: 
        thresh += 0.05 
        p1_x, p1_y, p2_x, p2_y, p_m = calc_perpendicular(a_x, b_x, m, b, thresh)
        data = get_positives_with_equation(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends, p1_x, p1_y, p_m)
        x_pt = data[2]
        y_pt = data[3]
        x.append(x_pt)
        y.append(y_pt)
    return x, y

def get_rates(orfs, genes, gene_ends):
    num_false_positives, num_true_positives = positive(genes, gene_ends)

    # TPR = TP/TP+FN
    # TP + FN = CDS
    tp_rate = num_true_positives*1.0/len(gene_ends)

    # FPR += FP/FP+TN
    # TN = orfs not meeting threshold or in cds
    fp_rate = num_false_positives*1.0/(len(orfs)-len(genes)-len(gene_ends)+num_true_positives)

    return (num_false_positives, num_true_positives, fp_rate, tp_rate)

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

# 3 ######################################################################

def len_thresh(orfs, gene_ends, threshold):
    genes = []
    for orf in orfs:
        seq, start, end = orf
        if len(seq) > threshold:
            genes.append(orf)

    return get_rates(orfs, genes, gene_ends)

def reach_len_thresh(orfs, gene_ends, threshold):
    num_false_positives, num_true_positives, fp_rate, tp_rate = len_thresh(orfs, gene_ends, threshold)
    tp_percent = tp_rate * 100
    print "LENGTH THRESHOLD of " + str(threshold) + " --> " + str(tp_percent) + "% true positive rate, " + str(num_true_positives) + " true positives, " + str(num_false_positives) + " false positives"
    print ""

# 4 ######################################################################

def mm_thresh(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs):
    genes = []
    for orf in orfs:
        seq, start, end = orf
        mm_score = markov_model.calc_mm_score(seq, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
        if mm_score > threshold:
            genes.append(orf)

    return get_rates(orfs, genes, gene_ends)

def reach_mm_thresh(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs):
    num_false_positives, num_true_positives, fp_rate, tp_rate = mm_thresh(orfs, gene_ends, threshold, k, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
    tp_percent = tp_rate * 100
    print "MM SCORE THRESHOLD of " + str(threshold) + " --> " + str(tp_percent) + "% true positive rate, " + str(num_true_positives) + " true positives, " + str(num_false_positives) + " false positives"
    print ""

# 5 ######################################################################

def plot_scatter(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends, thresh):
    plt.suptitle('Combining Length and MM Scores')
    plt.xlabel('ORF Length')
    plt.ylabel('Markov Model Score')

    trusted_gene_x, trusted_gene_y = plot_scatter_helper(k, trusted_genes, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)

    trusted_non_gene_x, trusted_non_gene_y = plot_scatter_helper(k, trusted_non_genes,trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)

    all_gene_x, all_gene_y, all_non_gene_x, all_non_gene_y = plot_scatter_distribute(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends)

    # plt.scatter(trusted_gene_x, trusted_gene_y, color='orange')
    # plt.scatter(trusted_non_gene_x, trusted_non_gene_y, color='blue')
    plt.scatter(all_gene_x, all_gene_y, color='orange')
    plt.scatter(all_non_gene_x, all_non_gene_y, color='blue')

    a_x, a_y, b_x, b_y, m, b = calc_median(trusted_gene_x, trusted_gene_y, trusted_non_gene_x, trusted_non_gene_y)

    p1_x, p1_y, p2_x, p2_y, p_m = calc_perpendicular(a_x, b_x, m, b, 0.24)
    num_false_positives, num_true_positives, fp_rate, tp_rate = get_positives_with_equation(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends, p1_x, p1_y, p_m)
    print "false positives:", num_false_positives, "true positives:", num_true_positives, "fpr:", fp_rate, "tpr:", tp_rate

    # median line 
    plt.annotate('A', (a_x, a_y), color='cyan')
    plt.annotate('B', (b_x, b_y), color='cyan')
    plt.plot([a_x, b_x], [a_y, b_y], color='cyan')

    # perpendicular line 
    plt.plot([p1_x, p2_x], [p1_y, p2_y], color='green')

    # draw lines x = 50, x = 1400
    plt.plot([50, 50], [-400, 300], color='red')
    plt.plot([1400, 1400], [-400, 300], color='magenta')

    plt.legend(['y = 0.0342799012441x - 0.63156280769 (median)', 'y - 11.6680657587 = m(x - 358.8) (perpendicular)', 'x = 50 (short ORF threshold)', 'x = 1400 (long ORF threshold)', 'genes', 'non-genes'], loc='lower right')
    
    plt.show()

def plot_scatter_helper(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs):
    x = []
    y = []
    for seq in orfs:
        x_pt = len(seq)
        y_pt = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
        x.append(x_pt)
        y.append(y_pt)
    return (x, y)

def plot_scatter_distribute(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends):
    gene_x = []
    gene_y = []
    non_gene_x = []
    non_gene_y = []

    for orf in orfs: 
        seq, start, end = orf 
        x_pt = len(seq)
        y_pt = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
        if end in gene_ends:
            gene_x.append(x_pt)
            gene_y.append(y_pt)
        else:
            non_gene_x.append(x_pt)
            non_gene_y.append(y_pt)
    
    return (gene_x, gene_y, non_gene_x, non_gene_y)

def calc_median(trusted_gene_x, trusted_gene_y, trusted_non_gene_x, trusted_non_gene_y):
    a_x = np.median(np.array(trusted_non_gene_x))
    a_y = np.median(np.array(trusted_non_gene_y))
    b_x = np.median(np.array(trusted_gene_x))
    b_y = np.median(np.array(trusted_gene_y))

    m = (a_y - b_y) / (a_x - b_x)
    b = a_y - m*a_x

    print "Median Line"
    print "A: (" + str(a_x) + ", " + str(a_y) + ")"
    print "B: (" + str(b_x) + ", " + str(b_y) + ")"
    print "slope: " + str(m)
    print "equation: y = " + str(m) + "x + " + str(b)
    print "" 

    return (a_x, a_y, b_x, b_y, m, b)

def calc_perpendicular(a_x, b_x, m, b, percent): 
    # intercept = % into slope
    pi_x = a_x + percent * (b_x - a_x)
    pi_y = m * pi_x + b
    
    # slope = negative recipricol
    p_m = -1 * np.reciprocal(m)

    # y - y1 = m(x - x1) --> y = mx -mx1 + y1
    # get 2 points on the line 
    p1_x = pi_x + 13
    p1_y = calc_y_coord(p1_x, pi_x, pi_y, p_m)
    p2_x = pi_x - 10
    p2_y = calc_y_coord(p2_x, pi_x, pi_y, p_m)

    print "Perpendicular Line"
    print "intercept: (" + str(pi_x) + ", " + str(pi_y) + ")"
    print "slope:", p_m
    print "equation: y -", pi_y, "= m(x -", pi_x, ")"
    print "point 1: (", p1_x, ",", p1_y, ")"
    print "point 2: (", p2_x, ",", p2_y, ")"
    print "" 

    return (p1_x, p1_y, p2_x, p2_y, p_m)

def calc_y_coord(x, x1, y1, m):
    return (m * x) - (m * x1) + y1

def get_positives_with_equation(k, orfs, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs, gene_ends, x1, y1, m):
    genes = []
    for orf in orfs: 
        seq, start, end = orf 
        x = len(seq)
        y = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes, trusted_genes_probs, trusted_non_genes_probs)
        y_thresh = calc_y_coord(x, x1, y1, m)
        if y > y_thresh:
            genes.append(orf)

    return get_rates(orfs, genes, gene_ends)

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

start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))