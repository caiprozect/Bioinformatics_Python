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

	def parse_refflat_line(self, sReadLine):
		sFlds = sReadLine.strip().split('\t')
		self.sRefID = sFlds[1]
		self.sGeneSymbol = sFlds[0]
		self.sChromID = sFlds[2][3:]
		self.sStrand = sFlds[3]
		self.nNumExons = int(sFlds[8])
		self.lnExStartPs = [int(startP) for startP in sFlds[9].strip(",").split(",")]
		self.lnExEndPs = [int(endP) for endP in sFlds[10].strip(",").split(",")]

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

def main(sFileName):
	listCRefSeq = parseWholeRefFlat(sFileName)
	print("Answer 1: {}".format(len(listCRefSeq)))
	listCRefSeq = parseNM(listCRefSeq)
	print("Answer 2: {}".format(len(listCRefSeq)))
	listCRefSeq = parseUniqueNM(listCRefSeq)
	print("Answer 3: {}".format(len(listCRefSeq)))
	listCRefSeq = sortByRefSeqID(listCRefSeq)

	listRefID = [cRefSeq.getRefID() for cRefSeq in listCRefSeq]
	listGeneSymbol = [cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]

	dfRefSeq = DataFrame({"Ref Seq ID": listRefID, "Gene Symbol": listGeneSymbol})
	dfRefSeq.to_excel("RefSeqID_GeneSymbol.xlsx", sheet_name="sheet1", index=False)

if __name__=="__main__":
	sInFile = "../data/refFlat.txt"
	main(sInFile)