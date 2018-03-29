import multiprocessing as mp #Python does not support proper threading, This is expensive to create than a thread
from collections import defaultdict #Auto initialization
from functools import reduce
from time import time
import os, re, psutil

def make_seq(infile, outfile): #Sanity check and Make one line seq file
	inFileH = open(infile, 'r')
	inFileH.readline() #Skip the first line
	inFileLines =  inFileH.read().upper().splitlines()
	inFileH.close()
	inFileNumbLine = len(inFileLines)
	filteredLines = [inLine.strip() for inLine in inFileLines if (bool(re.match('^[ACGTN]+$', inLine)))]
	assert(inFileNumbLine == len(filteredLines)), "\"{}\" file contains unknown nucleotide character".format(infile)
	outFileH = open(outfile, 'w')
	outFileH.write("".join(filteredLines))
	print("\t1. \"{}\" file passed sanity check".format(infile))
	print("\t2. Seq file \"{}\" has been created".format(outfile))
	print("\t=================================================")

def chunk_process(file, start, size, polyLens):
	lFreqs = []
	for polyLen in polyLens:
		partFreqDict = defaultdict(int)
		with open(file, "rb") as inFileH:
			inFileH.seek(start)
			sPartText = inFileH.read(size).upper()
			nPartTextLen = len(sPartText)
			for i in range(nPartTextLen-polyLen+1):
				partFreqDict[sPartText[i:i+polyLen]] += 1
		lFreqs.append(partFreqDict)
	#Endfor
	return lFreqs

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

def sum_dict_list(dictListA, dictListB):
	assert(len(dictListA) == len(dictListB)), "Two lists must have the same length"
	listLen = len(dictListA)
	summedList = []
	for i in range(listLen):
		summedList.append(sum_dicts(dictListA[i], dictListB[i]))
	return summedList

def freqForEach(file, polyLens):
	sFName = file
	sSFile = "../data/Seq.fa"
	make_seq(sFName, sSFile)

	nPolymerLen = polyLens #Counting monomers for this HW

	mpPool = mp.Pool() #Default option uses maximum number of cpu cores
	lJobs = []

	for ptrChunkStart, ptrChunkSize in make_chunk(sSFile):
		lJobs.append( mpPool.apply_async(chunk_process, (sSFile, ptrChunkStart, ptrChunkSize, nPolymerLen)) )

	lWFreqs = reduce((lambda x,y: sum_dict_list(x,y)), [job.get() for job in lJobs])

	#for sPolymer in lPolymers:
		#print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sPolymer.decode('utf-8'), wFreq[sPolymer], float(wFreq[sPolymer])/nWNumbPoly))

	mpPool.close()

	#print("\t=================================================")
	#print(psutil.Process().open_files())
	#assert(psutil.Process().open_files()==[]), "Bad: There are some open file handles!!!"
	#print("\tGood: All the file handles are properly closed!!")
	#print("\t---The End---")

	print("Frequency table for {} is generated\n\n".format(file))

	return lWFreqs

def main():
	polyLens = [1, 2]
	lChromNumb = [str(n) for n in range(1,23)] + ["X", "Y"]
	lInFiles = ["chr"+chromNumb+".fa" for chromNumb in lChromNumb]
	lInFilePaths = ["../data/chroms/"+inFileName for inFileName in lInFiles]
	lWFreqs = reduce((lambda x,y: sum_dict_list(x,y)), [freqForEach(each, polyLens) for each in lInFilePaths])
	monoFreq = lWFreqs[0]
	diFreq = lWFreqs[1]
	lMonomers = monoFreq.keys()
	lDimers = diFreq.keys()
	nWNumbMono = sum(monoFreq.values())
	nWNumbDi = sum(diFreq.values())

	print("\tMonomer stats\n")
	for sMonomer in lMonomers:
		print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sMonomer.decode('utf-8'), monoFreq[sMonomer], float(monoFreq[sMonomer])/nWNumbMono))
	print("\tDimer stats\n")
	for sDimer in lDimers:
		print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sDimer.decode('utf-8'), diFreq[sDimer], float(diFreq[sDimer])/nWNumbDi))

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("Took {} seconds to run".format(rtime))