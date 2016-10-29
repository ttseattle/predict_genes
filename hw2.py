import open_reading_frame 
import markov_model

MJANNASCHII_FILEPATH = "data/genome.fna"
ANNOTATIONS_FILEPATH = "data/annotation.gbff"
SHORT_THRESH = 50
LONG_THRESH = 1400
K = 3 

def main(): 
    # get data 
    seq = open_reading_frame.prep_mjannaschii(MJANNASCHII_FILEPATH) 
    gene_ends = open_reading_frame.prep_annotations(ANNOTATIONS_FILEPATH, seq)

    # homework problems 
    orfs = find_orfs(seq)
    short_orfs, long_orfs = find_short_and_long_orfs(orfs, SHORT_THRESH, LONG_THRESH)
    get_num_cdss(gene_ends)
    trusted_non_genes = remove_stop_codons(short_orfs)
    trusted_genes = remove_stop_codons(long_orfs)
    count_xy(K, trusted_genes, trusted_non_genes)
    print_summary(short_orfs, long_orfs, trusted_genes, trusted_non_genes, gene_ends)

# 1A #####################################################################

def find_orfs(seq): 
    print "NUMBER OF ORFS FOUND"
    orfs = []

    for offset in range(3):
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

def print_summary(short_orfs, long_orfs, trusted_genes, trusted_non_genes, gene_ends):
    print "SUMMARY FOR FIRST 5 ORFS OF LENGTH < " + str(SHORT_THRESH)
    print_summary_helper(short_orfs, trusted_genes, trusted_non_genes, gene_ends)
    print ""

    print "SUMMARY FOR FIRST 5 ORFS OF LENGTH > " + str(LONG_THRESH)
    print_summary_helper(long_orfs, trusted_genes, trusted_non_genes, gene_ends)
    print ""

def print_summary_helper(orfs, trusted_genes, trusted_non_genes, gene_ends):
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format("#", "Start", "Length", "MM Score", "CDS?")
    count = 1
    for orf in sorted(orfs, key=lambda tup: tup[1])[:5]:
        print_info(count, orf, trusted_genes, trusted_non_genes, gene_ends)
        count += 1

def print_info(count, orf, trusted_genes, trusted_non_genes, gene_ends): 
    seq, start, end = orf 
    mm_score = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes)
    cds = 'True' if end in gene_ends else 'False' 
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format(count, start, len(seq), mm_score, cds)

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