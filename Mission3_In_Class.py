from pandas import DataFrame
from collections import defaultdict
from functools import reduce
from time import time
import multiprocessing as mp 

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
		#Check
		self.nCDSStart = 0
		self.nCDSEnd = 0
		self.sCDSSeq = ""
		self.sORFCheck = ""

	def parse_refflat_line(self, sReadLine):
		sFlds = sReadLine.strip().split('\t')
		self.sRefID = sFlds[1]
		self.sGeneSymbol = sFlds[0]
		self.sChromID = sFlds[2][3:]
		self.sStrand = sFlds[3]
		self.nNumExons = int(sFlds[8])
		self.lnExStartPs = [int(startP) for startP in sFlds[9].strip(",").split(",")]
		self.lnExEndPs = [int(endP) for endP in sFlds[10].strip(",").split(",")]

		self.nCDSStart = int(sFlds[6])
		self.nCDSEnd = int(sFlds[7])

	def putExSeq(self, sExSeq):
		self.sExSeq = sExSeq
		self.nExSize = len(sExSeq)

	def putCDSSeq(self, sCDSSeq):
		self.sCDSSeq = sCDSSeq

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

	def putTrueCDS(self):
		if self.sStrand == '-':
			print("{}: is a - strand".format(self.sRefID))
			print(self.sCDSSeq[0:10])
			complNucl = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
			sComplCDSSeq = ""
			sCDSSeq = self.sCDSSeq
			for nucl in sCDSSeq:
				nucl = complNucl[nucl]
				sComplCDSSeq = nucl + sComplCDSSeq
			self.sCDSSeq = sComplCDSSeq
			print(self.sCDSSeq[-10:])

	def parseORF(self):
		#print("{}: finding ORF".format(self.sRefID))
		sExSeq = self.sExSeq
		sCDSSeq = self.sCDSSeq
		sStartC = 'ATG'
		lStopCs = ['TAG', 'TAA', 'TGA']
		
		if sCDSSeq[:3] == sStartC:
			if sCDSSeq[-3:] in lStopCs:
				if len(sCDSSeq)%3 == 0:
					lCodons = [sCDSSeq[i:i+3] for i in range(3, len(sCDSSeq)-3, 3)]
					if all((stopC not in lCodons) for stopC in lStopCs):
						nOffStart = sExSeq.find(sCDSSeq)
						assert(nOffStart >= 0), "Does not have matching CDS"
						nOffEnd = nOffStart + len(sCDSSeq)
						self.s5UTRSeq = sExSeq[:nOffStart]
						self.n5UTRSize = len(self.s5UTRSeq)
						self.sORFSeq = sExSeq[nOffStart:nOffEnd]
						self.nORFSize = len(self.sORFSeq)
						self.s3UTRSeq = sExSeq[nOffEnd:]
						self.n3UTRSize = len(self.s3UTRSeq)
						#print("{}: legitimate ORF".format(self.sRefID))
						assert(self.sORFSeq[:3] == 'ATG'), "Wrong Start Codon"
						assert(self.sORFSeq[-3:] in lStopCs), "Wrong End Codon"
					else:
						print("{}: internal stop codon".format(self.sRefID))
						self.sORFCheck = "INTERNAL STOP CODON"
				else:
					print("{}: not a multiple of 3".format(self.sRefID))
					self.sORFCheck = "No MULTIPLE OF 3"
			else:
				print("{}: does not end with STOP Codon".format(self.sRefID))
				self.sORFCheck = "NO STOP CODON"
		else:
			print("{}: does not start with ATG".format(self.sRefID))
			self.sORFCheck = "NO START CODON"

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

	def getCDSStart(self):
		return self.nCDSStart
	def getCDSEnd(self):
		return self.nCDSEnd
	def getCDSSeq(self):
		return self.sCDSSeq
	def getORFCheck(self):
		return self.sORFCheck
#End Class Definition

def genExSeqsPool(chromDict, listCRefSeq):
	sOrganism = list(chromDict.keys())[0] #Only deal with one organim
	lChromIDs = list(chromDict.values())[0] #Only deal with one organism

	mpPool = mp.Pool()
	lJobs = []

	for chID in lChromIDs:
		listMatch = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getChromID() == chID]
		lJobs.append( mpPool.apply_async(genExSeqs, (sOrganism, chID, listMatch)) )

	lProcRefSeq = reduce(lambda x, y: x+y, [job.get() for job in lJobs])

	mpPool.close()

	return lProcRefSeq


#Mission3 Functions

def genExSeqs(sOrganism, chID, listCRefSeq): #Does not account for strand +/-
	lProcRefSeq = []
	listMatch = listCRefSeq
	print("{}: Generating Exons".format(chID))
	sChFile = "../data/{}.chroms/chr{}.fa".format(sOrganism, chID)
	hChFile = open(sChFile, 'r')
	hChFile.readline() #Skip the first line
	sChromSeq = hChFile.read().upper().splitlines()
	hChFile.close()
	sChromSeq = "".join(sChromSeq)
	for letter in sChromSeq:
		assert(letter in ['A', 'C', 'G', 'T', 'N']), "{} chromosome File has unknown nucleotide character".format(chID)
	for cRefSeq in listMatch:
		nCDSStart = cRefSeq.getCDSStart()
		nCDSEnd = cRefSeq.getCDSEnd()
		print("{}: Trying to put exon".format(cRefSeq.getRefID()))
		sExSeq = ""
		lnExStartPs = cRefSeq.getExStartPs()
		lnExEndPs = cRefSeq.getExEndPs()
		nNumExons = cRefSeq.getNumExons()
		assert(len(lnExStartPs) == nNumExons), "{}: Number of exons does not match for starting positions".format(cRefSeq.getRefID())
		assert(len(lnExEndPs) == nNumExons), "{}: Number of exons does not match for end positions".format(cRefSeq.getRefID())
		lnExPs = list(zip(lnExStartPs, lnExEndPs))
		for ps in lnExPs:
			start = ps[0]
			end = ps[1]
			sExSeq = sExSeq + sChromSeq[start:end]
		cRefSeq.putExSeq(sExSeq)
		print("Exon: {}".format(cRefSeq.getExSeq()))
		print("{}: Exon successfully added".format(cRefSeq.getRefID()))

		print("{}: Trying to put CDS".format(cRefSeq.getRefID()))
		lnCDSPs = [ps for ps in lnExPs if ps[1] > nCDSStart]
		lnCDSPs = [ps for ps in lnCDSPs if ps[0] < nCDSEnd]
		lnCDSPs[0] = (nCDSStart, lnCDSPs[0][1])
		lnCDSPs[-1] = (lnCDSPs[-1][0] ,nCDSEnd)
		sCDSSeq = ""
		for ps in lnCDSPs:
			start = ps[0]
			end = ps[1]
			sCDSSeq = sCDSSeq + sChromSeq[start:end]
		cRefSeq.putCDSSeq(sCDSSeq)
		lProcRefSeq.append(cRefSeq)
	return lProcRefSeq

def takecareNegStrand(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.putTrueEx()
		cRefSeq.putTrueCDS()
		lProcCRefSeq.append(cRefSeq)
	return lProcCRefSeq

def genORFSeqs(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.parseORF()
		lProcCRefSeq.append(cRefSeq)
	return lProcCRefSeq

def onlyValidORF(listCRefSeq):
	sCheckF = "Mission3_CDS_Check.txt"
	hCheckF = open(sCheckF, 'w')
	noORFList = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getORFSize() == 0]
	for cRefSeq in noORFList:
		sRefID = cRefSeq.getRefID()
		sStrand = cRefSeq.getStrand()
		sORFCheck = cRefSeq.getORFCheck()
		sExSeq = cRefSeq.getExSeq()
		sCDSSeq = cRefSeq.getCDSSeq()
		print("{}\n{}\n\n{}\n\n{}\n\n{}\n\n\n\n\n\n".format(sRefID, sStrand, sORFCheck, sExSeq, sCDSSeq), file=hCheckF)
	hCheckF.close()
	return [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getORFSize() > 0]

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
	sOutFile = "Mission3_In_Class.txt"
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

	listRefID = [cRefSeq.getRefID() for cRefSeq in listCRefSeq]
	listGeneSymbol = [cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]
	list5UTRSize = [cRefSeq.get5UTRSize() for cRefSeq in listCRefSeq]
	listORFSize = [cRefSeq.getORFSize() for cRefSeq in listCRefSeq]
	list3UTRSize = [cRefSeq.get3UTRSize() for cRefSeq in listCRefSeq]

	dfRefSeq = DataFrame({"Ref Seq ID": listRefID, "Gene Symbol": listGeneSymbol,
		"5'UTR Size": list5UTRSize, "ORF Size": listORFSize, "3'UTR Size": list3UTRSize})
	dfRefSeq.to_excel("Mission3_In_Class.xlsx", sheet_name="sheet1", index=False)

	#Check Seqs
"""
>>>>>>> 9f7851d5b1238e43adec6a69567b13b2d7698efd
	checkFile = "Mission3_Seq_Check.txt"

	checkFileH = open(checkFile, "w")
	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		sExSeq = cRefSeq.getExSeq()
		sCDSSeq = cRefSeq.getCDSSeq()
		checkFileH.write(sRefID + "\n")
		checkFileH.write(sExSeq + "\n")
		checkFileH.write(sCDSSeq + "\n\n\n")

	checkFileH.close()

	checkNegStrand = "Mission3_Neg_Str_Chk.txt"
	checkNegH = open(checkNegStrand, "w")
	listNegs = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getStrand()=="-"]
	for cRefSeq in listNegs:
		sRefID = cRefSeq.getRefID()
		sStrand = cRefSeq.getStrand()
		sExSeq = cRefSeq.getExSeq()
		sCDSSeq = cRefSeq.getCDSSeq()
		checkNegH.write(sRefID + "\n")
		checkNegH.write(sStrand + "\n")
		checkNegH.write(sExSeq + "\n")
		checkNegH.write(sCDSSeq + "\n\n\n")

	checkNegH.close()
"""
"""
	sCheckFile = "Mission3_Check.txt"
	hCheckF = open(sCheckFile, 'w')
	shortORFs = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getORFSize() <= 12]
	for cRefSeq in shortORFs:
		sRefID = cRefSeq.getRefID()
		sStrand = cRefSeq.getStrand()
		sExSeq = cRefSeq.getExSeq()
		sCDSSeq = cRefSeq.getCDSSeq()
		print("{}\n{}\n{}\n{}".format(sRefID, sStrand, sExSeq, sCDSSeq), file=hCheckF)
	hCheckF.close()
"""
if __name__=="__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("{} seconds".format(rtime))