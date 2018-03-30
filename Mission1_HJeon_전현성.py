import multiprocessing as mp #Python does not support proper threading, This is expensive to create than a thread
from collections import defaultdict #Auto initialization
from functools import reduce
from time import time
import os, re, psutil

def numbToPtn(numb, polyLen):
	lNucls = ["A", "C", "G", "T"]
	nQ = numb
	nR = 0
	sPtn = ""
	for i in range(polyLen):
		nQ, nR = divmod(nQ, len(lNucls))
		sPtn = lNucls[nR] + sPtn
	return sPtn

def make_seq(infile): #Sanity check is done here
	inFileH = open(infile, 'r')
	inFileH.readline() #Skip the first line
	inFileLines =  inFileH.read().upper().splitlines()
	inFileH.close()
	inFileNumbLine = len(inFileLines)
	filteredLines = [inLine.strip() for inLine in inFileLines if (bool(re.match('^[ACGTN]+$', inLine)))]
	assert(inFileNumbLine == len(filteredLines)), "\"{}\" file contains unknown nucleotide character".format(infile)
	print("\t1. \"{}\" file passed sanity check".format(infile))
	print("\t2. Seq has been created")
	print("\t=================================================")
	return("".join(filteredLines))

def chunk_process(file, polyLen):
	lNucls = ["A", "C", "G", "T"]
	sSeqText = make_seq(file)
	tempMonoFreqs = defaultdict(int)
	tempPolyFreqs = defaultdict(int)
	monoFreqs = {}
	polyFreqs = {}
	expFreqs = {}
	windowSlideEnd = len(sSeqText) - polyLen + 1 
	for i in range(windowSlideEnd):
		if(i == windowSlideEnd):
			sPtn = sSeqText[i:i+polyLen]
			for letter in sPtn:
				tempMonoFreqs[letter] += 1
			tempPolyFreqs[sPtn] += 1
		else:
			sPtn = sSeqText[i:i+polyLen]
			tempMonoFreqs[sPtn[0]] += 1
			tempPolyFreqs[sPtn] += 1
	#Endfor
	for nucl in lNucls:
		monoFreqs[nucl] = tempMonoFreqs[nucl]
	for numb in range(4**polyLen):
		sPtn = numbToPtn(numb, polyLen)
		polyFreqs[sPtn] = tempPolyFreqs[sPtn]
		prob = 1
		for letter in sPtn:
			prob *= monoFreqs[letter]
		expFreqs[sPtn] = prob * (len(sSeqText)-polyLen) 
	return [monoFreqs, polyFreqs, expFreqs]

def sum_dicts(dictA, dictB):
	summedDict = defaultdict(int)
	unionKeys = set(dictA.keys()) | set(dictB.keys())
	for key in unionKeys:
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

def freqForEach(filelist, polyLen):
	mpPool = mp.Pool() #Default option uses maximum number of cpu cores
	lJobs = []

	for file in filelist:
		lJobs.append( mpPool.apply_async(chunk_process, (file, polyLen)) )

	lWFreqs = reduce((lambda x,y: sum_dict_list(x,y)), [job.get() for job in lJobs])

	#for sPolymer in lPolymers:
		#print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sPolymer.decode('utf-8'), wFreq[sPolymer], float(wFreq[sPolymer])/nWNumbPoly))

	mpPool.close()

	#print("\t=================================================")
	#print(psutil.Process().open_files())
	#assert(psutil.Process().open_files()==[]), "Bad: There are some open file handles!!!"
	#print("\tGood: All the file handles are properly closed!!")
	#print("\t---The End---")

	print("Frequency table is generated\n\n")

	return lWFreqs

def main():
	polyLen = 2 #Dealing with Dimers
	lChromNumb = [str(n) for n in range(1,23)] + ["X", "Y"]
	lInFiles = ["chr"+chromNumb+".fa" for chromNumb in lChromNumb]
	lInFilePaths = ["../data/chroms/"+inFileName for inFileName in lInFiles]
	lWFreqs = freqForEach(lInFilePaths, polyLen)
	monoFreq = lWFreqs[0]
	diFreq = lWFreqs[1]
	expFreq = lWFreqs[2]
	lMonomers = monoFreq.keys()
	lDimers = diFreq.keys()
	nWNumbMono = sum(monoFreq.values())
	nWNumbDi = sum(diFreq.values())
	nWexp = sum(expFreq.values())

	print("\tMonomer stats")
	for sMonomer in lMonomers:
		print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}".format(sMonomer, monoFreq[sMonomer], float(monoFreq[sMonomer])/nWNumbMono))
	print("\n\tDimer stats")
	for sDimer in lDimers:
		print("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f} \tFrequency of exp {0}: {3:f}".format(sDimer, diFreq[sDimer], float(diFreq[sDimer])/nWNumbDi, expFreq[sDimer]/nWexp))

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("Took {} seconds to run".format(rtime))