import multiprocessing as mp
import os
from collections import Counter
from functools import reduce

def chunk_process(start, size):
	c = Counter()
	with open(input_file) as f:
		f.seek(start)
		text = f.read(size)
		return Counter(text)

def make_chunk(file, size = 1024*1024*64):
	end = os.path.getsize(file)
	f = open(file, 'rb')
	cend = f.tell()
	while True:
		cstart = cend
		f.seek(size, 1)
		f.readline()
		cend = f.tell()
		yield cstart, cend - cstart
		if cend > end:
			f.close()
			break

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

pool = mp.Pool()
jobs = []

for cstart, csize in make_chunk(input_file):
	jobs.append( pool.apply_async(chunk_process, (cstart, csize)) )

result = reduce((lambda x,y: x+y), [job.get() for job in jobs])
print("Whole frequency dictionary: ", result)
print("# of A: ", result['A'])
print("# of C: ", result['C'])
print("# of G: ", result['G'])
print("# of T: ", result['T'])

pool.close()