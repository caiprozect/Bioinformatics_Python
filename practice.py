import multiprocessing as mp
import os
from collections import Counter
from functools import reduce

input_file = "human_genome"

def chunk_process(start, size):
	with open(input_file) as f:
		f.seek(start)
		text = f.read(size)
		return Counter(text)

def make_chunk(file, size = 1024*1024*1024):
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

pool = mp.Pool(4)
jobs = []

for cstart, csize in make_chunk(input_file):
	jobs.append( pool.apply_async(chunk_process, (cstart, csize)) )

result = reduce((lambda x,y: x+y), [job.get() for job in jobs])
print("# of A: ", result['A'])
print("# of C: ", result['C'])
print("# of G: ", result['G'])
print("# of T: ", result['T'])

pool.close()