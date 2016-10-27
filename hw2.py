import pdb

RNA_STOPS = ["UAG", "UAA", "UGA"]
DNA_STOPS = ["TAG", "TAA", "TGA"]
NUCLEOTIDES = ["A", "T", "C", "G"]
GENE_LEN = 1000

def findOrfs(seq): 
    genes, genome_ends = prep_annotations(seq)

    orfs = []
    for offset in range(3):
        orfs_at_offset = scan(seq, offset)
        orfs += orfs_at_offset
        print "number of ORFS at offset " + str(offset), len(orfs_at_offset)
    # orfs1 = scan(seq, 0)
    # orfs2 = scan(seq, 1)
    # orfs3 = scan(seq, 2)
    # orfs = orfs1+orfs2+orfs3

    print "total number of ORFS", len(orfs)

    count = 0 
    count_seq = 0

    orfs_l_50 = 0
    orfs_g_1400 = 0
    for orf in sorted(orfs, key=lambda tup: tup[1]): 
        seq = orf[0]
        start = orf[1]+1
        end = orf[2]+1

        if len(seq) < 50:
            orfs_l_50 += 1
        if len(seq) > 1400:
            orfs_g_1400 += 1

        # if len(seq) >= GENE_LEN:
        if end in genome_ends:
            # print "MATCHING " + str(start) + ".." + str(end) + " " + seq
            count = count + 1

        if seq in genes: 
            count_seq += 1
        
        # print str(start) + ".." + str(end) + " " + seq
    print "number of matching ORFS", count
    print "number of matching ORF sequences", count_seq

    print "less than 50", orfs_l_50
    print "greater than 1400", orfs_g_1400

    # for seq in genes:
    #     print seq, genes[seq]

def print_gene(seq, start, end):
    print str(start) + ".." + str(end) + " " + seq

def scan(seq, offset): 
    orfs = []

    orf_seq = ""
    idx = offset
    orf_start = idx 

    while idx <= len(seq): 
        triplet = seq[idx:idx + 3]
        orf_seq += triplet 

        if is_stop(triplet):
            orf_end = idx + 1 # want the index of the whole codon? 
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

                idx_range = content.split("..")
                start = int(idx_range[0])-1
                end = int(idx_range[1])

                gene = seq[start:end]
                genes[gene] = end 
                ends[end] = True
        
        # only process the first chromosome 
        if parts and parts[0] == "ORIGIN":
            break 
    
    print "number of annnotated genes", len(genes) 
    return (genes, ends)

def main(): 
    data = prep_mjannaschii() 
    # print data
    findOrfs(data)

main()