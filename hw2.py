import pdb

# define stop codons 
RNA_STOPS = ["UAG", "UAA", "UGA"]
DNA_STOPS = ["TAG", "TAA", "TGA"]
NUCS = ["A", "T", "C", 'G']
GENE_LEN = 1000

def findOrfs(seq): 
    genes, genome_ends = prep_annotations()

    orfs = []

    orfs1 = []
    orfs2 = []
    orfs3 = []

    # for offset in range(3):
    one = scan(orfs1, seq, 0)
    two = scan(orfs2, seq, 1)
    three = scan(orfs3, seq, 2)

    print "one", len(orfs1)
    print "two", len(orfs2)
    print "three", len(orfs3)
    
    print "number of ORFS", len(orfs)
    count = 0 
    for orf in sorted(orfs, key=lambda tup: tup[1]): 
        seq = orf[0]
        start = orf[1] + 1
        end = orf[2] + 1

        # if len(seq) >= GENE_LEN:
        if end in genome_ends:
            # print str(start) + ".." + str(end) + " " + seq
            count = count + 1
        # if seq in genes: 
        #     print seq 
    print "number of matching ORFS", count

def print_gene(seq, start, end):
    print str(start) + ".." + str(end) + " " + seq

def scan(orfs, seq, offset): 
    orf_start = offset
    orf = ""
    idx = orf_start

    count = 0 
    while idx <= len(seq): 
        # if count > 10: 
        #     break

        triplet = seq[idx:idx+3]
        orf = orf + triplet 
        if is_stop(triplet):
            # if len(orf)>=GENE_LEN:
                # print_gene(orf, orf_start, idx+3)
            orfs.append((orf, orf_start, idx+3))
            print orf
            count = count + 1
            orf = ""    #reset orf 
            # orf_start = idx + 3 
            # idx = orf_start
        idx = idx + 3 

# def scan(orfs, seq, offset):
#     orf_seq = ""
#     orf_start = offset
#     for i in range(offset, len(seq)+1, 3):  
#         triplet = seq[i:i+3]
#         orf_seq = orf_seq + triplet 
#         if is_stop(triplet):
#             orf = (orf_seq, orf_start, i+3)
#             orfs.append(orf) 

#             # print orf_seq
#             # break 

#             # reset orf
#             orf_seq = "" 
#             orf_start = i+3

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
        # pdb.set_trace()
        # only record the first long chromosome - stop at the beginning of the small ones 
        if line[0]==">":
            break

        for char in line: 
            if char == '\n':
                continue 
            elif char not in NUCS: 
                data += "T" 
            else:
                data += char

        # formatted_line = ""
        # for char in line: 
        #     # if char == ' ':
        #     #     continue 
        #     # elif char=='A' or char=='a':
        #     #     formatted_line = formatted_line + 'A'
        #     # elif char=='C' or char=='c':
        #     #     formatted_line = formatted_line + 'C'
        #     # elif char=='G' or char=='g':
        #     #     formatted_line = formatted_line + 'G'
        #     # else:    
        #     #     formatted_line = formatted_line + 'T'
        #     if char not in NUCS: 
        #         formatted_line += 'T'
        #     else:
        #         formatted_line += char 
        # data = data + formatted_line 
    
    return data 

def prep_annotations(): 
    seq = prep_mjannaschii() 
    file = open("data/annotation.gbff")

    fileContent = file.readlines() 

    ends = {}
    genes = {}

    for line in fileContent: 
        parts = line.split()
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
                start = int(idx_range[0])
                end = int(idx_range[1])
                gene = seq[start:end]
                genes[gene] = end 
                ends[end] = True
        # only process the first chromosome 
        elif parts and parts[0] == "ORIGIN":
            break 
    
    print "number of annnotated genes", len(genes) 
    return (genes, ends)

def main(): 
    data = prep_mjannaschii() 
    findOrfs(data)

    # data = "UUAAUGUGUCAUUGAUUAAG"
    # print findOrfs(data)

main()