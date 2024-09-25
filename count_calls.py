import sys

lines = {}
chr = ""

def printLines(lines):
    for line in lines:
       print(line+"\t"+str(lines[line]["CpG"])+"\t"+str(lines[line]["CHG"])+"\t"+str(lines[line]["CHH"])+"\t"+sys.argv[1])

for line in sys.stdin:
    split_line = line.split()
    gene_rec = "\t".join(split_line[:12])
    if split_line[0] != chr: #clear on new chrom to save memory and processing
        printLines(lines)
        lines = {}
        chr = split_line[0]
    if gene_rec not in lines:
        lines[gene_rec] = {"CpG": 0, "CHG": 0, "CHH": 0}
    context = split_line[18] 
    lines[gene_rec][context] += 1

printLines(lines)
