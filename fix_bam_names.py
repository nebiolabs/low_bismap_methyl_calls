import os
import os.path
import sys

def process_bam(bamname, aligner, ref):
    bamname = bamname.replace(".rg.deduplicated.sorted", "").replace(".md.sorted", "")
    name_for_params = bamname.rstrip(".bam").rsplit("_", 2)
    seqname = name_for_params[0]
    method = "unknown"
    mass = "unknown"
    if "SRR" in seqname:
        method = "bisulfite"
    elif "TS" in seqname:
        method = "bisulfite"
        mass = seqname.split("-")[0]
    elif "LB" in seqname:
        method = "emseq"
        mass = seqname.split("_")[0]
    newname = f"{seqname}.{mass}.{method}.{ref}.{aligner}.bam"
    return newname

def process_dir(dirname, dirpath):
    fields = dirname.split("_")
    ref = fields[1]
    aligner = fields[2]
    for bamname in os.listdir(dirpath):
        if bamname.endswith(".bam.bai"):
            continue
        newname = process_bam(bamname, aligner, ref)
        # assuming if X.bam exists, the index is X.bam.bai
        os.rename(os.path.join(dirpath, bamname), os.path.join(dirpath, newname))
        os.rename(os.path.join(dirpath, bamname+".bai"), os.path.join(dirpath, newname+".bai"))

bams_dir = sys.argv[1]

for dirname in os.listdir(bams_dir):
    process_dir(dirname, os.path.join(bams_dir, dirname))


