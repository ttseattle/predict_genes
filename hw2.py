import open_reading_frame 
import markov_model

MJANNASCHII_FILEPATH = "data/genome.fna"
ANNOTATIONS_FILEPATH = "data/annotation.gbff"
K = 3 

def main(): 
    # get data 
    seq = open_reading_frame.prep_mjannaschii(MJANNASCHII_FILEPATH) 
    gene_ends = open_reading_frame.prep_annotations(MJANNASCHII_FILEPATH, seq)

    # homework problems 
    orfs = find_orfs(seq)
    trusted_non_genes, trusted_genes = find_genes_and_non_genes(orfs, gene_ends)
    get_num_cdss(gene_ends)
    count_xy(K, trusted_non_genes)
    count_xy(K, trusted_genes)
    print_table(trusted_genes, trusted_non_genes)

# 1A #####################################################################

def find_orfs(seq): 
    orfs = []

    for offset in range(3):
        orfs_at_offset = open_reading_frame.scan(seq, offset)
        orfs += orfs_at_offset
        print "number of ORFS at offset " + str(offset), len(orfs_at_offset)

    print "total number of ORFS", len(orfs)
    print ""
    return orfs

# 1B & 1C ################################################################

def find_genes_and_non_genes(orfs, gene_ends):
    trusted_non_genes = []
    trusted_genes = []

    for orf in orfs: 
        seq, start, end = orf
        if len(seq) < 50:
            trusted_non_genes.append(orf)
        elif len(seq) > 1400:
            trusted_genes.append(orf)

    print "less than 50", len(trusted_non_genes)
    print "greater than 1400", len(trusted_genes)
    print ""
    return (trusted_non_genes, trusted_genes)

# 1D #####################################################################

def get_num_cdss(cdss):
    print "number of simple forward strand CDSs", len(cdss)
    print ""

# 1E #####################################################################

def count_xy(k, training):
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

def print_table(trusted_genes, trusted_non_genes):
    print "First 5 ORFS of length less than 50"
    print_table_helper(trusted_non_genes, trusted_genes, trusted_non_genes)
    print ""

    print "First 5 ORFS of length greater than 1400"
    print_table_helper(trusted_genes, trusted_genes, trusted_non_genes)
    print ""

def print_table_helper(seqs, trusted_genes, trusted_non_genes):
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format("#", "Start", "Length", "Markov Model", "CDS?")
    count = 1
    for orf in sorted(seqs, key=lambda tup: tup[1])[:5]:
        print_info(count, orf, trusted_genes, trusted_non_genes)
        count += 1

def print_info(count, orf, trusted_genes, trusted_non_genes): 
    seq, start, end = orf 
    mm_score = markov_model.calc_mm_score(seq, K, trusted_genes, trusted_non_genes)
    cds = 'True' if seq in trusted_genes else 'False' 
    
    print "{:<5} {:<10} {:<10} {:<18} {:<10}".format(count, start, len(seq), mm_score, cds)

# ERROR CHECKING #########################################################

def check_genes(orfs, gene_ends):
    count = 0
    for orf in orfs: 
        subseq, start, end = orf
        if end in gene_ends:
            count += 1
    print "matching ends", count 
    print ""

# EXECUTE PROGRAM ########################################################
main()