import xlwt

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
		sFlds = sReadLine.split('\t')
		self.sRefID = sFlds[1]
		self.sGeneSymbol = sFlds[0]
		self.sChromID = sFlds[2][3:]
		self.sStrand = sFlds[3]
		self.nNumExons = int(sFlds[8])
		self.lnExStartPs = [int(startP) for startP in sFlds[9].strip(",").split(",")]
		self.lnExEndPs = [int(endP) for endP in sFlds[10].strip(",").split(",")]

	def getRefID(self):
		return sRefID
	def getGeneSymbol(self):
		return sGeneSymbol
	def getChromID(self):
		return sChromID
	def getStrand(self):
		return sStrand
	def getNumExons(self):
		return nNumExons
	def getExStartPs(self):
		return lnExStartPs
	def getExEndPs(self):
		return lnExEndPs

def parseWholeRefFlat(sFileName):
	wholeCRefSeq = []
	fileH = open(sFileName, 'r')
	for sLine in fileH:
		cRefSeq = RefSeq()
		cRefSeq.parse_refflat_line(sLine)
		wholeCRefSeq.append(cRefSeq)
	return wholeCRefSeq

def parseNM(listCRefSeq):
	chromList = [str(n) for n in range(1,23)] + ['X', 'Y']
	filteredList = []
	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		sHead = sRefID.split("_")[0]
		sChromID = cRefSeq.getChromID()
		if (sHead == "NM") and (sChromID in chromList):
			filteredList.append(cRefSeq)
	return filteredList

def parseUniqueNM(sortedCRefSeqs):
	uniqueList = []
	nSize = len(sortedCRefSeqs)
	nIdx = 0
	while nIdx < nSize:
		sCurr = sortedCRefSeqs[nIdx].getRefID()
		sNext = sortedCRefSeqs[nIdx+1].getRefID()
		if (sCurr == sNext):
			nIdx += 2
			while sortedCRefSeqs[nIdx].getRefID() == sCurr:
				nIdx += 1
		

def sortByRefSeqID(listCRefSeq):
	listCRefSeq.sort(key = (lambda cRefSeq: int(cRefSeq.getRefID().split("_")[-1])))
	return listCRefSeq

def main(sFileName):
	listCRefSeq = parseWholeRefFlat(sFileName)
	print("Answer 1: {}".format(len(listCRefSeq)))
	listCRefSeq = parseNM(listCRefSeq)
	print("Answer 2: {}".format(len(listCRefSeq)))
	listCRefSeq = sortByRefSeqID(listCRefSeq)
	listCRefSeq = parseUniqueNM(listCRefSeq)
	print("Answer 3: {}".format(len(listCRefSeq)))

	xslH = xlwt.Workbook(encoding='utf-8')
	xslSheet = xslH.add_sheet("RefID, Gene Symbol")

	xslSheet.write(0, 0, "RefSeqID")
	xslSheet.write(0, 1, "Gene Symbol")

	rowNum  = 1

	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		sGeneSymbol = cRefSeq.getGeneSymbol()
		xslSheet.write(rowNum, 0, sRefID)
		xslSheet.write(rowNum, 1, sGeneSymbol)
		rowNum += 1

	xslH.save("RefSeqID_GeneSymbol")

if __name__=="__main__":
	sInFile = "./data/refFlat.txt"
	main(sInFile)