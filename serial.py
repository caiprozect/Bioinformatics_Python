from collections import Counter

fname = "./chroms/chr1.fa"

with open(fname, 'r') as f:
	f.readline()
	text = f.read().upper()
	count = Counter(text)
	whole = count['A']+count['C']+count['G']+count['T']
	print("A: ", count['A'], " freq A: ", count['A']/whole)
	print("T: ", count['T'], " freq T: ", count['T']/whole)
	print("C: ", count['C'], " freq C: ", count['C']/whole)
	print("G: ", count['G'], " freq G: ", count['G']/whole)