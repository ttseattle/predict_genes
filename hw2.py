import pdb
import math 
import markovmodel

RNA_STOPS = ["UAG", "UAA", "UGA"]
DNA_STOPS = ["TAG", "TAA", "TGA"]
NUCLEOTIDES = ["A", "C", "G", "T"]
GENE_LEN = 1000
K = 3 

def findOrfs(): 
    seq = prep_mjannaschii() 

    # find all ORFs within the sequence 
    orfs = []
    for offset in range(3):
        orfs_at_offset = scan(seq, offset)
        orfs += orfs_at_offset
        print "number of ORFS at offset " + str(offset), len(orfs_at_offset)
    print "total number of ORFS", len(orfs)

    # get the real genes from the annotations
    gene_ends, genes = prep_annotations(seq)

    count = 0
    count2 = 0
    for orf in orfs: 
        # pdb.set_trace()
        subseq, start, end = orf
        if end in gene_ends:
            # pdb.set_trace()
            count += 1
            if subseq in genes: 
                count2 += 1
            # else: 
            #     pdb.set_trace()
    print "matching ends", count 
    print "matching genes", count2

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

    markovmodel.part_e(K, trusted_genes)
    markovmodel.part_e(K, trusted_non_genes)
    
    print_table(trusted_genes, trusted_non_genes)

def print_table(trusted_genes, trusted_non_genes):
    for orf in sorted(trusted_non_genes, key=lambda tup: tup[1])[:5]:
        seq, start, end = orf
        print_info(seq, start, end, trusted_genes, trusted_non_genes)

    for orf in sorted(trusted_genes, key=lambda tup: tup[1])[:5]:
        seq, start, end = orf
        print_info(seq, start, end, trusted_genes, trusted_non_genes)

def print_info(seq, start, end, trusted_genes, trusted_non_genes): 
    print ""
    print seq 
    print "Start Index", start 
    print "End Index", end 
    print "Length", len(seq)
    markovmodel.calc_mm_score(seq, K, trusted_genes, trusted_non_genes)
    print "CDS?", "T" if seq in trusted_genes else "F"
    print "-------------------"
    print ""

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

def prep_mjannaschii(): 
    file = open("data/genome.fna")
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

def prep_annotations(seq): 
    file = open("data/annotation.gbff")
    fileContent = file.readlines() 

    ends = {}
    genes = {}
    for line in fileContent: 
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
                gene = seq[start:end+1]

                if end==10004:
                    print gene

                genes[gene] = True 
        
        # only process the first chromosome 
        if parts and parts[0] == "ORIGIN":
            break 

    return (ends, genes)

def main(): 
    # data = prep_mjannaschii() 
    # findOrfs(data)
    findOrfs()

    # training = ["AAAAGCTA"]
    # print markovmodel.condit_prob("A","AA",2,training)
    # print markovmodel.condit_prob("G","AA",2,training)
    # print markovmodel.condit_prob("A","CT",2,training)
    # print markovmodel.uncondit_prob("AA",2,training)
    # print markovmodel.uncondit_prob("AG",2,training)    

main()