import sys
import argparse

# Parse argument
parser = argparse.ArgumentParser()
parser.add_argument("-r","--rmgff", help="Reformatted .gff from repeatmasker (by reformatRMgff.py)", type=str)
parser.add_argument("-s","--size", help="Short TEs size for all TEs class except simple repeat, low complexity and satellite", type=int)
parser.add_argument("-f","--mfile", help="Mapfile (.tsv)",type=str)
args = parser.parse_args()

if len(sys.argv) <= 1:
	parser.print_help()
	sys.exit()

if (args.rmgff is None):
	parser.print_help()
	sys.exit("Reformatted repeatmasker gff is missing")

if (args.size is not None) and (args.mfile is not None):
	parser.print_help()
	sys.exit("Please provide either size or mapfile.")

gff = args.rmgff

if args.mfile is None:
	s = 100
else:
	s = None

mapfile = args.mfile


# if mapfile.tsv is not available, use defalut map
if s is not None:
	m = {
		"Unknown":s,
		"SINE":s,
		"SINE?":s,
		"LINE":s,
		"LTR" :s,
		"DNA" :s,
		"RC"  :s,
		"rRNA":s,
		"Simple_repeat":0,
		"Low_complexity":0,
		"Satellite":0,
		"snRNA":0
	}
else:
	m = {}
	with open(mapfile,"r") as f:
		for line in f:
			c,t = line.rstrip().split("\t")
			m[c] = int(t)
	# Assume SINE? = SINE
	m["SINE?"] = m["SINE"]

# Find how many line to skip in gff (faster for long file)
cnt = 0
with open(gff,"r") as f:
	for line in f:
			cnt +=1
			if not line.startswith("#"):
				cnt -=1
				break

# Check if size of TEs samller than lower threshold
with open(gff,"r") as f:
	for i in range(cnt):
		next(f)
	for line in f:
		col = line.rstrip().split("\t")
		generalClass = col[2].split("/")[0]
		teSize = int(col[4]) - int(col[3])
		try:
			m[generalClass]
		except KeyError:
			m[generalClass] = 0
			sys.stderr.write("'" + generalClass + "' is not in the mapfile. Skip filtering " + generalClass)
		if teSize < m[generalClass]:
			col[8] =  col[8] + ";shortTE=T"
		else:
			col[8] = col[8] + ";shortTE=F"
		print(*col,sep = "\t")




