from pandas import DataFrame
from collections import defaultdict

class RefSeq:
	def __init__(self):
		self.sRefID = ""
		self.sGeneSymbol = ""
		self.sChromID = ""
		self.sStrand = ""
		self.nNumExons = 0
		self.lnExStartPs = []
		self.lnExEndPs = []
		#Mission3 Intermediate
		self.sExSeq = "" #It is a concatenated + strand DNA seq!
		self.sORFSeq = ""
		#Mission3
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

	def parseORF(self):
		

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
#End Class Definition

#Mission3 Functions
def genExSeqs(chromDict, listCRefSeq): #Does not accunt for strand +/-
	sOrganism = list(chromDict.keys())[0] #Only deal with one organim
	lChromIDs = list(chromDict.values())[0] #Only deal with one organism
	dictCRefByChromID = {}
	for chID in lChromIDs:
		dictCRefByChromID[chID] = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getChromID() == chID]
	for chID in lChromIDs:
		print("{}: Generating Exons".format(chID))
		sChFile = "../data/{}.chroms/chr{}.fa".format(sOrganism, chID)
		hChFile = open(sChFile, 'r')
		hChFile.readline() #Skip the first line
		sChromSeq = hChFile.read().upper().splitlines()
		hChFile.close()
		sChromSeq = "".join(sChromSeq)
		for letter in sChromSeq:
			assert(letter in ['A', 'C', 'G', 'T', 'N']), "{} chromosome File has unknown nucleotide character".format(chID)
		for cRefSeq in dictCRefByChromID[chID]:
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
				print("{} Exon: {} successfully added".format(cRefSeq.getRefID(), sExSeq))

def genORFSeqs(listCRefSeq):
	[cRefSeq.parseORF() for cRefSeq in listCRefSeq]

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

	return uniqueList

def sortByRefSeqID(listCRefSeq):
	listCRefSeq.sort(key = (lambda cRefSeq: int(cRefSeq.getRefID().split("_")[-1])))
	return listCRefSeq

def main():
	sFileName = "../data/refFlat.txt"
	listCRefSeq = parseWholeRefFlat(sFileName)
	print("Answer 1: {}".format(len(listCRefSeq)))
	listCRefSeq = parseNM(listCRefSeq)
	print("Answer 2: {}".format(len(listCRefSeq)))
	listCRefSeq = parseUniqueNM(listCRefSeq)
	print("Answer 3: {}".format(len(listCRefSeq)))
	listCRefSeq = sortByRefSeqID(listCRefSeq)

	chromDict = {'hg38': ([str(i) for i in range(1,23)]+['X','Y'])}
	genExSeqs(chromDict, listCRefSeq)
	genORFSeqs(listCRefSeq)	

if __name__=="__main__":
	main()