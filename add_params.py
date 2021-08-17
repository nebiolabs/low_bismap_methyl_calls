#!/usr/bin/env python
import sys

masses = {}

def getMass(line, index):
    l = line.split()
    if len(l) < index+1:
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

with open(sys.argv[1]) as f:
    for line in f:
        out.write(line+"\t"+sys.argv[4]+"\t"+sys.argv[5]+"\t"+sys.argv[6]+"\t"+sys.argv[7]+"\t"+getMass(line, sys.argv[8])+"\n")