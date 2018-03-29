infile = "human_genome"
outfile = "sequence"

f = open(outfile, 'w')
f.close()

with open(infile, 'r') as f:
	with open(outfile, 'a') as out:
		for line in f:
			if line[0] != '>':
				out.write(line)