import sys
import re

if len(sys.argv) != 3:
    sys.exit("Usage: rmGffAddClass2Attr.py <repeatmasker.gff> <repeatmasker.out> \
             \nReformt repeatmasker gff. Type column -> General repeat class; Attribute -> Target_start:(pos); Target_end:(pos); ID:(repeat-family)")

rmgff = sys.argv[1]
rmout = sys.argv[2]


# Find out the repeat:repeatClass pair

classD = {}
with open(rmout,"r") as f:
    for i in range(3): # Skip .out header
        next(f)
    for line in f:
        [_, _, _, _, _, _, _, _, _, repeat, repeatClass, _, _, _, _] = line.rstrip().split()[0:15]
        classD[repeat] = repeatClass

# Rewrite the attr in repeatmasker gff
print("##gff-version 3")
with open(rmgff, "r") as f:
    for line in f:
        if line.startswith("#"): # skip header
            continue
        else:
            [seqid, source, T, start, end, score, strand, phase, remark] = line.rstrip().split("\t")
            family = re.findall("Motif:(.*)\"", remark)[0]
            s,e = remark.split()[2:]
            c = classD[family]
            nremark = "Tstart="+ s + ";Tend=" + e +";ID=" + family
            print(*[seqid,source,c,start,end,score,strand,phase,nremark],sep="\t")

