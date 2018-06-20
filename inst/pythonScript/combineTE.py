import sys
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)


# arguments
if len(sys.argv) <= 1:
    sys.exit("Usage: combineTE.py [reformatted .gff] [gap size, default=150]")

gapSize = 150
if len(sys.argv) == 2:
    gff = sys.argv[1]
else:
    gff = sys.argv[1]
    gapSize = int(sys.argv[2])

# test input
# gff = "/home/user/wongwaiyee/scratch/hydra/cleaning_pipeline/sample/hsym.fa.reformat.gff"

# Main
# Initial parameters dict
P = {
    "pcol": [],
    "pEnd": 0,
    "pFamily": "",
    "pCstart": 0,
    "pCend": 0,
    "plabel": ""
}

D = nested_dict()


def update_pcol(c, label=""):
    attr = c[8].split(";")
    attrD = {}
    for i in attr:
        k,v = i.split("=")
        attrD[k] = v

    P["pcol"] = c
    P["pEnd"] = int(c[4])
    P["pFamily"] = attrD["ID"]
    P["pCstart"] = int(attrD["Tstart"])
    P["pCend"] = int(attrD["Tend"])
    P["plabel"] = label


print("##gff-version 3")
with open(gff, "r") as f:
    for line in f:
        # If current TE is far away from the last TE
        col = line.rstrip().split("\t")
        if int(col[3])-P["pEnd"] > gapSize:
            if P["pcol"]:  # Make sure not the first line
                print(*P["pcol"], sep="\t")
            update_pcol(c=col, label="")
        else:
            # Extract attribute
            cattrD = {}
            cattr = col[8].split(";")
            for i in cattr:
                k, v = i.split("=")
                cattrD[k] = v

            # If current TE close to last TE, but not the same family
            if P["pFamily"] != cattrD["ID"]:
                if P["pcol"]:  # Make sure not the first line
                    print(*P["pcol"], sep="\t")
                update_pcol(c=col, label="")
            # If current TE close to last TE and belong to same family
            # Since has compared the family, pcol must not be empty in the following block
            else:
                o = min(int(P["pCend"]), int(cattrD["Tend"])) - max(int(P["pCstart"]), int(cattrD["Tstart"]))
                if o > 0:  # Consensus position overlap
                    print(*P["pcol"], sep="\t")
                    update_pcol(c=col, label="")
                else:
                    # Group the two TE by adding group attr
                    # If the previous TE is already in a group
                    if P["plabel"]:
                        print(*P["pcol"], sep="\t")
                        update_pcol(c=col, label=P["plabel"])  # keep using the same label
                        P["pcol"][8] = P["pcol"][8] + ";" + grouplabel # Add attribute
                        # Previous TE is the first element in the group, add the attribute before print
                    # Previous TE has no label yet
                    else:
                        # Check label
                        if D[col[0]][col[2]][cattrD["ID"]]:  # Use a new label for the TE pair
                            D[col[0]][col[2]][cattrD["ID"]] += 1
                            groupID = D[col[0]][col[2]][cattrD["ID"]]
                           # D[col[0]][col[2]][cattrD["ID"]] = groupID
                        else:  # If this is the first family in this chrom
                            D[col[0]][col[2]][cattrD["ID"]] = 1
                            groupID = 1
                        grouplabel = "TEgroup=" + col[0] + "|" + col[2] + "|" + cattrD["ID"] + "|" + str(groupID)
                        P["pcol"][8] = P["pcol"][8] + ";" + grouplabel
                        print(*P["pcol"], sep="\t")
                        update_pcol(c=col, label=grouplabel)
                        P["pcol"][8] = P["pcol"][8] + ";" + grouplabel  # Add attribute


# Print the last row
print(*P["pcol"], sep="\t")
