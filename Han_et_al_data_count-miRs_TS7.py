# miRNA counting for Han et al. data
# Peter Y. Wang 2023
# Bartel Lab, Whitehead Institute/MIT

import sys
from collections import defaultdict as dd
import re

# Build look up
speciesRE = re.compile("([a-z]{3})\\-")
miRdict = {}
with open("./miR_Family_Info_TS7.txt", "r") as d:
    headline = d.readline()
    for l in d:
        items = l.strip("\n").split("\t")
        miR = items[3]
        species = speciesRE.match(miR).groups()[0]
        if species == "hsa":
            miRdict[items[4][:19]] = miR

# Get counts
cts = dd(list)
unmapped = []
totaln = 0
with open(sys.argv[1], "r") as r:
    for l in r:
        totaln += 1

        readseq = l.strip("\n").replace("T","U")
        barcode = readseq[:4] + readseq[-4:]
        insert = readseq[4:-4]
        insert19 = insert[:19]

        if len(insert) >= 20 and len(insert) <= 30 and insert19 in miRdict.keys():
            cts[insert19].append(barcode)
        else:
            unmapped.append(readseq)

# Output
print("total reads # = " + str(totaln))
print("mapped reads # = " + str(totaln-len(unmapped)) + " (" + str((1-len(unmapped)/totaln)*100) + "%)")

with open(sys.argv[2]+".csv", "w") as c, open(sys.argv[2]+"_unmapped.txt", "w") as u:
    c.write("miR,seq19,counts\n")
    for ent in cts.keys():
        c.write(",".join([
                miRdict[ent],
                ent,
                str(len(cts[ent]))  # no barcode collapsing
                ]) + "\n")
    for s in unmapped:
        u.write(s+"\n")
