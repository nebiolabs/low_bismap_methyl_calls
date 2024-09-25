#!/usr/bin/env python
import sys

for line in sys.stdin:
    rec = line.split("\t")[int(sys.argv[2])].strip()
    rec = rec.split(".")
    print(line.strip()+"\t"+rec[1]+"\t"+rec[2]+"\t"+rec[3]+"\t"+rec[4]+"\t"+sys.argv[1]+"\t"+rec[0])
