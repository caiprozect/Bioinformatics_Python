import numpy as np
from functools import reduce

def ptnToNumb(sPtn):
	lNucls = ['A', 'C', 'G', 'T']
	nNumb = 0
	nPtnLen = len(sPtn)
	for i in range(nPtnLen):
		nNumb += 4**(nPtnLen-1-i) * lNucls.index(sPtn[i])
	return nNumb

def numbToPtn(nNumb, nK):
	lNucls = ['A', 'C', 'G', 'T']
	sPtn = ""
	nCurrQ = nNumb
	nCurrR = 0
	for i in range(nK):
		nCurrR = nCurrQ % 4
		sPtn = lNucls[nCurrR] + sPtn
		nCurrQ = int(nCurrQ/4)
	return sPtn

def cptFreq(file, nK):
	f = open(file, 'r')
	f.readline()
	lines = f.read().upper().splitlines()
	f.close()
	sText = "".join(lines)
	nparrFreqArray = np.zeros(4**nK)
	nTextLen = len(sText)
	for i in range(nTextLen - nK + 1):
		sPtn = sText[i:i+nK]
		try:
			nNumb = ptnToNumb(sPtn)
		except:
			continue
		nparrFreqArray[nNumb] += 1
	print("{} done for this file".format(file))
	return nparrFreqArray

def processOrganChromDict(organChromDict, nK, outfile):
	with open(outfile, 'w') as f:	
		seqNameList = organChromDict.keys()
		for seqName in seqNameList:
			chromNumbs = organChromDict[seqName]
			chromFileList = ["../data/{}.chroms/chr{}.fa".format(seqName, numb) for numb in chromNumbs]
			freqs = np.array([cptFreq(file, nK).astype(np.int64) for file in chromFileList])
			freqs = np.sum(freqs, axis=0)
			flds = "{:16}{:>16}{:>16}{:>16}\n"
			entries = "{:16}{:16d}\n"
			f.write("{} chromosomes statistics\n".format(seqName))
			f.write(flds.format("K-mer", "Counts", "Frequencies", "Expectations"))
			for i in range(freqs.size):
				kmer = numbToPtn(i, nK)
				f.write(entries.format(kmer, freqs[i]))
			f.write("\n")

def main():
	#dictChromCat = {'hg38': ([str(i) for i in range(1,23)]+['X','Y']), 'galGal3': ([str(i) for i in range(1,29)]+['32','W','Z']),
	#				'dm3': ['2R', '2L', '3R', '3L', '4', 'X'], 'ce10': ['I', 'II', 'III', 'IV', 'V', 'X']}
	dictChromCat = {'dm3': ['2R', '2L', '3R', '3L', '4', 'X']}
	K = 2
	outFile = "Mission1_Verify.txt"
	processOrganChromDict(dictChromCat, K, outFile)

if __name__=="__main__":
	main()