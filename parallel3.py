#Need sanity check!!!

import multiprocessing as mp
import numpy as np
from collections import defaultdict, Counter
from functools import reduce
import os

def chunk_process(file, start, size, length):
	partFreqDict = defaultdict(int)
	with open(file, "rb") as f:
		f.seek(start)
		sPartText = f.read(size).upper()
		nPartTextLen = len(sPartText)
		for i in range(nPartTextLen-length+1):
			partFreqDict[sPartText[i:i+length]] += 1
		return Counter(partFreqDict)

def make_chunk(file, size = 1024*1024*64):
	fileEnd = os.path.getsize(file)
	f = open(file, 'rb')
	f.readline()
	chunkEnd = f.tell()
	while True:
		chunkStart = chunkEnd
		f.seek(size, 1)
		f.readline()
		chunkEnd = f.tell()
		yield chunkStart, chunkEnd - chunkStart
		if chunkEnd > fileEnd:
			f.close()
			break

def main():
	sFName = "./chroms/chr1.fa"
	nPolymerLen = 1

	pool = mp.Pool()
	jobs = []

	for ptrChunkStart, ptrChunkSize in make_chunk(sFName):
		jobs.append( pool.apply_async(chunk_process, (sFName, ptrChunkStart, ptrChunkSize, nPolymerLen)) )

	wFreq = reduce((lambda x,y: x+y), [job.get() for job in jobs])
	nWNumbPoly = sum(wFreq.values()) - wFreq['\n'] - wFreq['N']
	lPolymers = wFreq.keys()
	lPolymers = [polymer for polymer in lPolymers if '\n'.encode('utf-8') not in polymer]
	lPolymers = [polymer for polymer in lPolymers if 'N'.encode('utf-8') not in polymer]

	for sPolymer in lPolymers:
		print("Number of "+sPolymer.decode('utf-8')+": ", wFreq[sPolymer], "\tFrequency of "+sPolymer.decode('utf-8')+": ", float(wFreq[sPolymer])/nWNumbPoly)

	pool.close()

main()