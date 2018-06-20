import sys

if len(sys.argv) != 3:
	sys.exit("Wrong argument. getSeq.py <.fasta> <seqID>")

targetid = str(sys.argv[2])

# Flag
seq2print = False

with open(sys.argv[1], "r") as f:
	for line in f:
		if not seq2print:
			if line.startswith(">"):
				#print(line.lstrip(">"))
				if line.rstrip().lstrip(">") == targetid:
					print(line.rstrip())
					seq2print = True
					continue
				else:
					continue
			else:
				continue
		else:  # seq2print == Ture
			if not line.startswith(">"):
				print(line.rstrip())
			else:
				break
