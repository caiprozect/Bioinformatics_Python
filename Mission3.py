from pandas import DataFrame
from collections import defaultdict
import multiprocessing as mp 
import re
from functools import reduce
import numpy as np

class RefSeq:
	def __init__(self):
		self.sRefID = ""
		self.sGeneSymbol = ""
		self.sChromID = ""
		self.sStrand = ""
		self.nNumExons = 0
		self.lnExStartPs = []
		self.lnExEndPs = []
		#Mission3
		self.sExSeq = "" #It is a concatenated + strand DNA seq!
		self.sORFSeq = ""
		self.s5UTRSeq = ""
		self.s3UTRSeq = ""
		self.nExSize = 0
		self.n5UTRSize = 0
		self.nORFSize = 0
		self.n3UTRSize = 0

	def parse_refflat_line(self, sReadLine):
		sFlds = sReadLine.strip().split('\t')
		self.sRefID = sFlds[1]
		self.sGeneSymbol = sFlds[0]
		self.sChromID = sFlds[2][3:]
		self.sStrand = sFlds[3]
		self.nNumExons = int(sFlds[8])
		self.lnExStartPs = [int(startP) for startP in sFlds[9].strip(",").split(",")]
		self.lnExEndPs = [int(endP) for endP in sFlds[10].strip(",").split(",")]

	def putExSeq(self, sExSeq):
		self.sExSeq = sExSeq
		self.nExSize = len(sExSeq)

	def putTrueEx(self):
		if self.sStrand == '-':
			print("{}: is a - strand".format(self.sRefID))
			print(self.sExSeq[0:10])
			complNucl = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
			sComplExSeq = ""
			sExSeq = self.sExSeq
			for nucl in sExSeq:
				nucl = complNucl[nucl]
				sComplExSeq = nucl + sComplExSeq
			self.sExSeq = sComplExSeq
			print(self.sExSeq[-10:])

	def parseORF(self):
		print("{}: finding ORF".format(self.sRefID))
		sExSeq = self.sExSeq
		nIdxStart = sExSeq.find('ATG')
		nCheck = 0
		if nIdxStart >= 0:
			lStopCs = ['TAG', 'TAA', 'TGA']
			nIdxStop = nIdxStart + 3
			while nIdxStop < (self.nExSize - 3 + 1):
				sCurrCodon = sExSeq[nIdxStop:nIdxStop+3]
				if sCurrCodon in lStopCs: 
					print("{}: ORF found".format(self.sRefID))
					nIdxStop = nIdxStop + 3
					self.s5UTRSeq = sExSeq[0:nIdxStart]
					self.n5UTRSize = len(self.s5UTRSeq)
					self.sORFSeq = sExSeq[nIdxStart:nIdxStop]
					self.nORFSize = len(self.sORFSeq)
					self.s3UTRSeq = sExSeq[nIdxStop:]
					self.n3UTRSize = len(self.s3UTRSeq)
					print(self.sORFSeq)
					nCheck = 1
					break
				nIdxStop += 3
			if nCheck == 0:
				print("{}: no valid stop codon".format(self.sRefID))
				print(self.sStrand)
				print(self.sChromID)
				print(self.lnExStartPs)
				print(self.lnExEndPs)
				print(self.sExSeq)
		else:
			print("{}: no valid start codon".format(self.sRefID))

	def getRefID(self):
		return self.sRefID
	def getGeneSymbol(self):
		return self.sGeneSymbol
	def getChromID(self):
		return self.sChromID
	def getStrand(self):
		return self.sStrand
	def getNumExons(self):
		return self.nNumExons
	def getExStartPs(self):
		return self.lnExStartPs
	def getExEndPs(self):
		return self.lnExEndPs
	def getExSeq(self):
		return self.sExSeq
	def getORFSeq(self):
		return self.sPORFSeq
	def getExSize(self):
		return self.nExSize
	def get5UTRSize(self):
		return self.n5UTRSize
	def getORFSize(self):
		return self.nORFSize
	def get3UTRSize(self):
		return self.n3UTRSize
	def get5UTRSeq(self):
		return self.s5UTRSeq
	def get3UTRSeq(self):
		return self.s3UTRSeq
#End Class Definition

#Mission3 Functions
def genExSeqsPool(chromDict, listCRefSeq):
	mpPool = mp.Pool()
	lJobs = []

	sOrganism = list(chromDict.keys())[0] #Only deal with one organim
	lChromIDs = list(chromDict.values())[0] #Only deal with one organism
	dictByChID = {}
	for chID in lChromIDs:
		dictByChID[chID] = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getChromID() == chID]

	for chID in lChromIDs:
		lJobs.append( mpPool.apply_async(genExSeqs, (sOrganism, chID, dictByChID[chID])) )

	lProcCRefSeq = reduce((lambda x,y: x + y), [job.get() for job in lJobs])

	mpPool.close()

	return lProcCRefSeq

def genExSeqs(sOrganism, chID, listCRefSeq): #Does not account for strand +/-
	lProcRefSeq = []
	print("{}: Generating Exons".format(chID))
	sChFile = "../data/{}.chroms/chr{}.fa".format(sOrganism, chID)
	hChFile = open(sChFile, 'r')
	hChFile.readline() #Skip the first line
	sChromSeq = hChFile.read().upper().splitlines()
	hChFile.close()
	sChromSeq = "".join(sChromSeq)
	for letter in sChromSeq:
		assert(letter in ['A', 'C', 'G', 'T', 'N']), "{} chromosome File has unknown nucleotide character".format(chID)
	for cRefSeq in listCRefSeq:
		print("{}: Trying to put exon".format(cRefSeq.getRefID()))
		sExSeq = ""
		lnExStartPs = cRefSeq.getExStartPs()
		lnExEndPs = cRefSeq.getExEndPs()
		nNumExons = cRefSeq.getNumExons()
		assert(len(lnExStartPs) == nNumExons), "{}: Number of exons does not match for starting positions".format(cRefSeq.getRefID())
		assert(len(lnExEndPs) == nNumExons), "{}: Number of exons does not match for end positions".format(cRefSeq.getRefID())
		lnExPs = zip(lnExStartPs, lnExEndPs)
		for ps in lnExPs:
			start = ps[0]
			end = ps[1]
			sExSeq = sExSeq + sChromSeq[start:end]
		cRefSeq.putExSeq(sExSeq)
		lProcRefSeq.append(cRefSeq)
		print("Exon: {}".format(cRefSeq.getExSeq()))
		print("{}: Exon successfully added".format(cRefSeq.getRefID()))
	return lProcRefSeq

def takecareNegStrand(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.putTrueEx()
		lProcCRefSeq.append(cRefSeq)
	return lProcCRefSeq

def genORFSeqs(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.parseORF()
		lProcCRefSeq.append(cRefSeq)
	return lProcCRefSeq

def onlyValidORF(listCRefSeq):
	return [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getORFSize() != 0]

def isoformFiltering(listCRefSeq):
	setSeenGS = set([])
	listFiltered = []
	for cRefSeq in listCRefSeq:
		sGS = cRefSeq.getGeneSymbol()
		if sGS not in setSeenGS:
			listFiltered.append(cRefSeq)
			setSeenGS.add(sGS)
	return listFiltered

#Mission2 Functions
def parseWholeRefFlat(sFileName):
	wholeCRefSeq = []
	fileH = open(sFileName, 'r')
	for sLine in fileH:
		cRefSeq = RefSeq()
		cRefSeq.parse_refflat_line(sLine)
		wholeCRefSeq.append(cRefSeq)
	fileH.close()
	return wholeCRefSeq

def parseNM(listCRefSeq):
	chromList = [str(n) for n in range(1,23)] + ['X', 'Y']
	excludedRefID = set([])
	excludedChromID = set([])
	filteredList = []
	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		sHead = sRefID.split("_")[0]
		sChromID = cRefSeq.getChromID()
		if (sHead == "NM") and (sChromID in chromList):
			filteredList.append(cRefSeq)
		else:
			if (sHead != "NM"):
				excludedRefID.add(sRefID)
			else:	
				excludedChromID.add(sChromID)
	outH = open("Excluded", 'w')
	print("not NM_", file=outH)
	[print(sRefID, file=outH) for sRefID in excludedRefID]
	print("\n\nnot Chrom Num", file=outH)
	[print(sChromID, file=outH) for sChromID in excludedChromID]
	outH.close()
	return filteredList

def parseUniqueNM(listCRefSeq):
	refDict = defaultdict(list)
	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		refDict[sRefID].append(cRefSeq)

	uniqueList = [cRefSeq for cRefSeq in listCRefSeq if len(refDict[cRefSeq.getRefID()]) == 1]
	overlapList = set([cRefSeq.getRefID() for cRefSeq in listCRefSeq if len(refDict[cRefSeq.getRefID()]) != 1])

	outH = open("NotUnique", 'w')
	print("not Unique", file=outH)
	[print(sRefID, file=outH) for sRefID in overlapList]
	outH.close()

	return uniqueList

def sortByRefSeqID(listCRefSeq):
	listCRefSeq.sort(key = (lambda cRefSeq: int(cRefSeq.getRefID().split("_")[-1])))
	return listCRefSeq

def main():
	sFileName = "../data/refFlat.txt"
	sOutFile = "Mission3.txt"
	hOutF = open(sOutFile, 'w')
	listCRefSeq = parseWholeRefFlat(sFileName)
	print("Answer 1: {}".format(len(listCRefSeq)), file=hOutF)
	listCRefSeq = parseNM(listCRefSeq)
	print("Answer 2: {}".format(len(listCRefSeq)), file=hOutF)
	listCRefSeq = parseUniqueNM(listCRefSeq)
	print("Answer 3: {}".format(len(listCRefSeq)), file=hOutF)

	chromDict = {'hg38': ([str(i) for i in range(1,23)]+['X','Y'])}
	listCRefSeq = genExSeqsPool(chromDict, listCRefSeq)
	listCRefSeq = takecareNegStrand(listCRefSeq)
	listCRefSeq = genORFSeqs(listCRefSeq)

	listCRefSeq = onlyValidORF(listCRefSeq)
	print("Answer 4: {}".format(len(listCRefSeq)), file=hOutF)

	listCRefSeq = sortByRefSeqID(listCRefSeq)

	listCRefSeq = isoformFiltering(listCRefSeq)
	print("Answer 5: {}".format(len(listCRefSeq)), file=hOutF)

	hOutF.close()	

if __name__=="__main__":
	main()