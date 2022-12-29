#!/usr/bin/env python
import sys

masses = {}

def getMass(line, index):
    l = line.split()
    if index is None or len(l) < index+1:
        return "unknown_mass"
    else:
        nameStr = l[index]
        if nameStr in masses:
            return masses[nameStr]
        else:
            return "unknown_mass"

with open(sys.argv[3]) as f:
    for line in f:
        mass = line.split(",")
        masses[mass[0]] = mass[1]

out = open(sys.argv[2], "w")

index = None

if len(sys.argv) >= 9:
    try:
        index = int(sys.argv[8])
    except ValueError: #argument is non-numeric
        pass # leaving it as None

with open(sys.argv[1]) as f:
    for line in f:
        out.write(line+"\t"+sys.argv[4]+"\t"+sys.argv[5]+"\t"+sys.argv[6]+"\t"+sys.argv[7]+"\t"+getMass(line, index)+"\n")