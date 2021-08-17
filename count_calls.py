import sys

lines = {}
out = open(sys.argv[1], "w")
chr = ""


def printLines(lines):
    for line in lines:
        out.write(line+"\t"+str(lines[line])+"\t"+sys.argv[3]+"\n")

with open(sys.argv[2]) as f:
    for line in f:
        if line.split()[0] != chr: #clear on new chrom to save memory and processing
            printLines(lines)
            lines = {}
            chr = line.split()[0]
        if line.strip() not in lines:
            lines[line.strip()] = 0
        lines[line.strip()] += 1
        #print(line.strip())

    printLines(lines)
