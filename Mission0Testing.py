import multiprocessing as mp 
from collections import defaultdict #Auto initialization
from functools import reduce
from time import time
import os, re, psutil

def make_seq(infile, outfile): #Sanity check and Make one line seq file
	inFileH = open(infile, 'r')
	inFileH.readline() #Skip the first line
	inFileLines =  inFileH.read().upper().splitlines()
	inFileH.close()
	outFileH = open(outfile, 'w')
	outFileH.write("".join(inFileLines))
	outFileH.close()
	print("\tSeq file has been created")
	print("\t=================================================")

def chunk_process(file, start, size, polyLen):
	partFreqDict = defaultdict(int)
	with open(file, "rb") as inFileH:
		inFileH.seek(start)
		sPartText = inFileH.read(size).upper()
		nPartTextLen = len(sPartText)
		for i in range(nPartTextLen-polyLen+1):
			sPtn = sPartText[i:i+polyLen].decode('utf-8')
			if 'N' not in sPtn:
				partFreqDict[sPtn] += 1
		#Endfor
	return partFreqDict

def make_chunk(file, polyLen, size = 1024*1024*64):
	fileEnd = os.path.getsize(file)
	inFileH = open(file, 'rb')
	chunkEnd = inFileH.tell()
	while True:
		chunkStart = (chunkEnd - polyLen + 1) if (chunkEnd - polyLen + 1 >= 0) else chunkEnd
		inFileH.seek(chunkStart)
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
		summedDict[key] = dictA[key] + dictB[key]
	#EndFor
	return summedDict

def main():
	sFName = "../data/chroms/chr1.fa"
	sSFile = "../data/Seq.fa"
	make_seq(sFName, sSFile)

	nPolymerLen = 1 #Counting monomers for this HW

	mpPool = mp.Pool() #Default option uses maximum number of cpu cores
	lJobs = []

	for ptrChunkStart, ptrChunkSize in make_chunk(sSFile, nPolymerLen):
		lJobs.append( mpPool.apply_async(chunk_process, (sSFile, ptrChunkStart, ptrChunkSize, nPolymerLen)) )

	wFreq = reduce((lambda x,y: sum_dicts(x,y)), [job.get() for job in lJobs])
	lPolymers = wFreq.keys()
	nWNumbPoly = sum(wFreq.values())

	for sPolymer in lPolymers:
		print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sPolymer, wFreq[sPolymer], float(wFreq[sPolymer])/nWNumbPoly))

	mpPool.close()

	print("\t=================================================")
	#print(psutil.Process().open_files())
	#assert(psutil.Process().open_files()==[]), "Bad: There are some open file handles!!!"
	#print("\tGood: All the file handles are properly closed!!")
	print("\t---The End---")

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("Took {} seconds to run".format(rtime))