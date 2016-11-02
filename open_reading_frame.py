import pdb
import math 
import markov_model

DNA_STOPS = ["TAG", "TAA", "TGA"]
NUCLEOTIDES = ["A", "C", "G", "T"]
GENE_LEN = 1000

def scan(seq, offset): 
    orfs = []

    orf_seq = ""
    idx = offset
    orf_start = idx 
    while idx <= len(seq): 
        triplet = seq[idx:idx + 3]
        orf_seq += triplet 
        if is_stop(triplet):
            orf_end = idx + 3
            orfs.append((orf_seq, orf_start, orf_end))
            #reset orf 
            orf_seq = ""    
            orf_start = idx + 3 
        idx += 3

    return orfs 

def is_stop(triplet):
    for dna_stop in DNA_STOPS: 
        if dna_stop==triplet:
            return True
    return False 

def prep_mjannaschii(filepath): 
    file = open(filepath)
    fileContent = file.readlines()[1:]  # skip first line

    data = ""
    for line in fileContent: 
        # only record the first long chromosome - stop at the beginning of the small ones 
        if line[0]==">":
            break

        for char in line: 
            if char == '\n' or char == ' ':
                continue 
            elif char not in NUCLEOTIDES: 
                data += "T" 
            else:
                data += char
    
    return data 

def prep_annotations(filepath, seq): 
    file = open(filepath)
    file_content = file.readlines() 

    ends = {}
    for line in file_content: 
        parts = line.split()

        # gene indices are located on lines labeled with CDS
        if parts and parts[0] == "CDS":
            content = parts[1]

            # ignore complements 
            if "complement" in content:
                continue 
            elif ".." in content:
                # for now, ignore < / >
                if "<" in content or ">" in content: 
                    continue
                # record at what index the sequence ended 
                idx_range = content.split("..")
                start = int(idx_range[0])
                end = int(idx_range[1])
                ends[end] = True
        
        # only process the first chromosome 
        if parts and parts[0] == "ORIGIN":
            break 

    return ends

# false-positive: not a gene, but says it is 
# true-positive: a gene, and says is a gene
def positive(orfs, gene_ends):
    num_false_positives = 0
    num_true_postives = 0
    for orf in orfs:
        seq, start, end = orf 
        if end in gene_ends:
            num_true_postives += 1
        else: 
            num_false_positives +=1
    return (num_false_positives, num_true_postives)

def len_thresh(orfs, gene_ends, threshold):
    genes = []
    for orf in orfs:
        seq, start, end = orf 
        if len(seq) > threshold:
            genes.append(orf)
    
    num_false_positives, num_true_positives = positive(genes, gene_ends)
    print "number false positives", num_false_positives
    print "number true positives", num_true_positives
    ratio = num_true_positives/(len(genes)*1.0)
    print "true positive rate", ratio

# MJANNASCHII_FILEPATH = "data/genome.fna"
# ANNOTATIONS_FILEPATH = "data/annotation.gbff"

# seq = prep_mjannaschii(MJANNASCHII_FILEPATH) 
# gene_ends = prep_annotations(ANNOTATIONS_FILEPATH, seq)

# # homework problems 
# orfs = hw2.find_orfs(seq)
# len_thresh(orfs, gene_ends, 1000)
