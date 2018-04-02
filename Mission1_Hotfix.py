import multiprocessing as mp 
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
	#Endfor
	return sPtn

def make_seq(infile): #Sanity check is done here
	inFileH = open(infile, 'rb')
	inFileH.readline() #Skip the first line
	inFileLines =  inFileH.read().upper().splitlines()
	inFileH.close()
	#Let's do sanity check in other lines since this takes up too much time
	#inFileNumbLine = len(inFileLines)
	#filteredLines = [inLine.strip() for inLine in inFileLines if (bool(re.match('^[ACGTN]+$', inLine)))]
	#assert(inFileNumbLine == len(filteredLines)), "\"{}\" file contains unknown nucleotide character".format(infile)
	print("\t{} Seq has been created".format(infile))
	print("\t=================================================")
	return (("".encode('utf-8')).join(inFileLines)).decode('utf-8')

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

	#Let's do sanity check here
	assert(set(tempMonoFreqs.keys()).issubset(set(["A", "C", "G", "T", "N"]))), "{} file has unknown nucleotide character {}".format(file, tempMonoFreqs.keys())

	for nucl in lNucls:
		monoFreqs[nucl] = tempMonoFreqs[nucl]
	#Endfor
	for numb in range(4**polyLen):
		sPtn = numbToPtn(numb, polyLen)
		polyFreqs[sPtn] = tempPolyFreqs[sPtn]
	#Endfor
	return [monoFreqs, polyFreqs]

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
	sOutFile = "Mission1_Hotfix.txt"
	polyLen = 2 #Dealing with Dimers
	dictChromCat = {'hg38': ([str(i) for i in range(1,23)]+['X','Y']), 'galGal3': ([str(i) for i in range(1,29)]+['32','W','Z']),
					'dm3': ['2R', '2L', '3R', '3L', '4', 'X'], 'ce10': ['I', 'II', 'III', 'IV', 'V', 'X']}

	outFileH = open(sOutFile, 'w')
	for chromVer, chromNumb in dictChromCat.items():
		outFileH.write("{} Statistics\n".format(chromVer))
		lChromNumb = chromNumb
		lInFiles = ["chr"+chromNumb+".fa" for chromNumb in lChromNumb]
		lInFilePaths = ["../data/{}.chroms/".format(chromVer)+inFileName for inFileName in lInFiles]
		lWFreqs = freqForEach(lInFilePaths, polyLen)
		monoFreq = lWFreqs[0]
		diFreq = lWFreqs[1]
		lMonomers = monoFreq.keys()
		lDimers = diFreq.keys()
		nWNumbMono = sum(monoFreq.values())
		nWNumbDi = sum(diFreq.values())
		
		expFreq = {}
		for key in diFreq.keys():
			expFreq[key] = 1
			for letter in key:
				expFreq[key] *= float(monoFreq[letter])/nWNumbMono
		#Endfor

		outFileH.write("\tMonomer stats\n")
		for sMonomer in lMonomers:
			outFileH.write("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f}\n".format(sMonomer, monoFreq[sMonomer], float(monoFreq[sMonomer])/nWNumbMono))
		#Endfor
		outFileH.write("\n\tDimer stats\n")
		for sDimer in lDimers:
			outFileH.write("\tNumber of {0}: {1:d} \tFrequency of {0}: {2:f} \tFrequency of exp {0}: {3:f}\n".format(sDimer, diFreq[sDimer], float(diFreq[sDimer])/nWNumbDi, expFreq[sDimer]))
		#Endfor
		outFileH.write("\n\n")
	#Endfor
	outFileH.close()

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("Took {} seconds to run".format(rtime))