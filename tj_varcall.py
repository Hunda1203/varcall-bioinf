import re
from collections import defaultdict
import csv

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

def parse_cigar(cigar):
    return [(int(length), op) for length, op in CIGAR.findall(cigar)]

read_type = input("Input read type (ecoli or lambda): ")
path_ref = read_type + ".fasta"
path_sam = read_type + ".sam"
min_map_quality = 20
min_reads = 3
min_share = 0.5


if __name__ == "__main__":
    ref = ''.join(open(path_ref).read().splitlines()[1:])
    subs = defaultdict(lambda: defaultdict(int))
    inserts = defaultdict(lambda: defaultdict(int))
    dels = defaultdict(int)
    
    total_reads = [0] * len(ref)
    
    with open(path_sam, "r") as path:
        for line in path:
            if line[0] == "@":
                continue
            fields = line.rstrip("\n").split("\t")
            
            flag = int(fields[1])
            if flag & 4 or flag & 256 or flag & 2048:
                continue
            
            map_quality = int(fields[4])
            if map_quality < min_map_quality:
                continue
            
            pos = int(fields[3]) - 1
            
            cigar = fields[5]
            if cigar == "*":
                continue
            
            
            seq = fields[9]
            if seq == "*":
                continue
            ref_pos = pos
            read_pos = 0
            
            cigar = parse_cigar(cigar)
            
            
            for i in range(len(cigar)):
                if cigar[i][1] == "=":
                    for j in range(cigar[i][0]):
                        total_reads[ref_pos + j] += 1
                    ref_pos += cigar[i][0]
                    read_pos += cigar[i][0]
                    
                elif cigar[i][1] == "X":
                    for j in range(cigar[i][0]):
                        total_reads[ref_pos + j] += 1
                        subs[ref_pos + j][seq[read_pos+j]] += 1
                    ref_pos += cigar[i][0]
                    read_pos += cigar[i][0]
                    
                elif cigar[i][1] == "I":
                    insert = seq[read_pos:read_pos+cigar[i][0]]
                    inserts[ref_pos][insert] += 1
                    read_pos += cigar[i][0]
                    
                elif cigar[i][1] == "D":
                    for j in range(cigar[i][0]):
                        total_reads[ref_pos + j] += 1
                        dels[ref_pos + j] += 1
                    ref_pos += cigar[i][0]
                    
                elif cigar[i][1] == "S":
                    read_pos += cigar[i][0]
                    
        calls = []
        
        
        for pos in subs:
            mut = list(subs[pos].keys())[0]
        
            for base in subs[pos]:
                if subs[pos][base] > subs[pos][mut]:
                    mut = base
        
            if subs[pos][mut] >= min_reads and total_reads[pos] > 0 and subs[pos][mut] / total_reads[pos] >= min_share:
                calls.append(("X", pos, mut))
        
        
        for pos in inserts:
            inserted = list(inserts[pos].keys())[0]
        
            for seq in inserts[pos]:
                if inserts[pos][seq] > inserts[pos][inserted]:
                    inserted = seq
        
            if pos < len(total_reads):
                depth = total_reads[pos]
            else:
                depth = total_reads[pos - 1]
        
            if inserts[pos][inserted] >= min_reads and depth > 0 and inserts[pos][inserted] / depth >= min_share:
                calls.append(("I", pos, inserted))
        
        
        for pos in dels:
            if dels[pos] >= min_reads and total_reads[pos] > 0 and dels[pos] / total_reads[pos] >= min_share:
                calls.append(("D", pos, "-"))
    
    
        calls.sort(key=lambda x: x[1])

        with open(read_type+"_calls.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerows(calls)
            