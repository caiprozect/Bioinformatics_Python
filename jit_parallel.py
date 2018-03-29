import numba
import os
import numpy as np
from collections import Counter

@numba.jit(nopython=True, parallel=True)
def chunk_process(start, size):
	c = Counter()
	with open(input_file) as f:
		f.seek(start)
		text = f.read(size)
		c = Counter(text)
		return np.array([c['A'], c['C'], c['G'], c['T']])

@numba.jit(nopython=True, parallel=True)
def make_chunk(file, size = 1024*1024*64):
	positions = []
	end = os.path.getsize(file)
	f = open(file, 'rb')
	cend = f.tell()
	while True:
		cstart = cend
		f.seek(size, 1)
		f.readline()
		cend = f.tell()
		positions.append((cstart, cend-cstart))
		if cend > end:
			return positions

infile = "human_genome"
outfile = "sequence"
input_file = outfile

f = open(outfile, 'w')
f.close()

with open(infile, 'r') as f:
	with open(outfile, 'a') as out:
		for line in f:
			if line[0] != '>':
				out.write(line)

positions = make_chunk(input_file)

@numba.jit(nopython=True, parallel=True)
def frequency_table():
	result = np.array([0,0,0,0])
	for cstart, csize in positions:
		result += chunk_process(cstart, cend)
	print("Whole frequency dictionary: ", result)
	print("# of A: ", result['A'])
	print("# of C: ", result['C'])
	print("# of G: ", result['G'])
	print("# of T: ", result['T'])

frequency_table()