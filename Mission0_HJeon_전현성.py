import multiprocessing as mp #Python does not support proper threading, This is expensive to create than a thread
from collections import defaultdict #Auto initialization
from functools import reduce
import os, re, psutil

def make_seq(infile, outfile): #Sanity check and Make one line seq file
	inFileH = open(infile, 'r')
	inFileH.readline() #Skip the first line
	open(outfile, 'w') #Empty the file if it already exists
	outFileH = open(outfile, 'a')
	inFileLines =  inFileH.read().upper().splitlines()
	inFileH.close()
	for inLine in inFileLines:
		inLine = inLine.strip()
		assert(bool(re.match('^["A","C","G","T","N"]+$', inLine))), "File contains unknown nucleotide character"
		outFileH.write(inLine)
	#EndFor
	outFileH.close()
	print("\t1. File passed sanity check")
	print("\t2. Seq file has been created")
	print("\t=================================================")

def chunk_process(file, start, size, length):
	partFreqDict = defaultdict(int)
	with open(file, "rb") as inFileH:
		inFileH.seek(start)
		sPartText = inFileH.read(size).upper()
		nPartTextLen = len(sPartText)
		for i in range(nPartTextLen-length+1):
			partFreqDict[sPartText[i:i+length]] += 1
	return partFreqDict

def make_chunk(file, size = 1024*1024*64):
	fileEnd = os.path.getsize(file)
	inFileH = open(file, 'rb')
	chunkEnd = inFileH.tell()
	while True:
		chunkStart = chunkEnd
		inFileH.seek(size, 1)
		chunkEnd = inFileH.tell()
		yield chunkStart, chunkEnd - chunkStart
		if chunkEnd > fileEnd:
			inFileH.close()
			break
	#EndWhile

def sum_dicts(dictA, dictB):
	summedDict = defaultdict(int)
	unionKeys = set(dictA.keys()) | set(dictB.keys())
	for key in unionKeys:
		if 'N'.encode('utf-8') not in key:
			summedDict[key] = dictA[key] + dictB[key]
	#EndFor
	return summedDict

def main():
	sFName = "./chroms/chr1.fa"
	sSFile = "Seq.fa"
	make_seq(sFName, sSFile)

	nPolymerLen = 1 #Counting monomers for this HW

	mpPool = mp.Pool() #Default option uses maximum number of cpu cores
	lJobs = []

	for ptrChunkStart, ptrChunkSize in make_chunk(sSFile):
		lJobs.append( mpPool.apply_async(chunk_process, (sSFile, ptrChunkStart, ptrChunkSize, nPolymerLen)) )

	wFreq = reduce((lambda x,y: sum_dicts(x,y)), [job.get() for job in lJobs])
	lPolymers = wFreq.keys()
	nWNumbPoly = sum(wFreq.values())

	for sPolymer in lPolymers:
		print("\tNumber of "+sPolymer.decode('utf-8')+": ", wFreq[sPolymer], "\tFrequency of "+sPolymer.decode('utf-8')+": ", float(wFreq[sPolymer])/nWNumbPoly)

	mpPool.close()

	print("\t=================================================")
	assert(psutil.Process().open_files()==[]), "Bad: There are some open file handles!!!"
	print("\tGood: All the file handles are properly closed!!")
	print("\t---The End---")

main()