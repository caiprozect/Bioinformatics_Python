from time import time
from scipy import stats
from pandas import DataFrame
import numpy as np
from collections import defaultdict
from functools import reduce
import multiprocessing as mp 
import pickle

class miRNA:
	def __init__(self):
		self.sSeq = ""
		self.sName = ""

	def parse_info(self, sName, sSeq):
		self.sSeq = sSeq
		self.sName = sName

	def getName(self):
		return self.sName
	def getSeq(self):
		return self.sSeq

#Mission 5 functions
def fill_miRNA(sMiRNAFile):
	hMiFile = open(sMiRNAFile, 'r')
	sMiData = hMiFile.read().split('>')
	hMiFile.close()

	listFilledCMiRNA = []

	for item in sMiData:
		flds = item.strip().split('\n')
		if len(flds) == 2:
			sSpecy = flds[0].split('-')[0]
			if sSpecy == "hsa":
				cMiRNA = miRNA()
				sName = flds[0].strip().split()[-1]
				sSeq = flds[1].strip().upper()
				sRTSeq = reverse_transcript(sSeq)
				cMiRNA.parse_info(sName, sRTSeq)
				listFilledCMiRNA.append(cMiRNA)

	return listFilledCMiRNA

def reverse_transcript(sSeq):
	dictRT = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}

	sRTSeq = ""
	for letter in sSeq[1:]:
		rtLetter = dictRT[letter]
		sRTSeq = rtLetter + sRTSeq

	return sRTSeq

def bonferroni_correction(listCMotif):
	listCorrCMotif = []
	bonferroniCoeff = len(listCMotif)

	for cMotif in listCMotif:
		fPValue = cMotif.getPValue()
		fCorrPValue = min(1., fPValue * bonferroniCoeff)
		cMotif.putPValue(fCorrPValue)
		listCorrCMotif.append(cMotif)

	return listCorrCMotif

def fill_target_info(listCMotif, listCMiRNA):
	dictCMotif = {}

	for cMotif in listCMotif:
		dictCMotif[cMotif.getMotif()] = cMotif

	for cMiRNA in listCMiRNA:
		sName = cMiRNA.getName()
		sSeq = cMiRNA.getSeq()
		s7mer_m8_seq = sSeq
		s7mer_A1_seq = sSeq[1:] + "A"

		if s7mer_m8_seq in dictCMotif:
			cM8Motif = dictCMotif[s7mer_m8_seq]
			cM8Motif.put_miRNA_Info(sName, "7mer-m8")
			dictCMotif[s7mer_m8_seq] = cM8Motif

		if s7mer_A1_seq in dictCMotif:
			cA1Motif = dictCMotif[s7mer_A1_seq]
			cA1Motif.put_miRNA_Info(sName, "7mer-A1")
			dictCMotif[s7mer_A1_seq] = cA1Motif

	return list(dictCMotif.values())


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
		#NMD
		self.bNMD = False

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

	def checkNMD(self):
		nAfterStop = self.nExSize - self.n3UTRSize
		if self.sStrand == "+":
			nLastEx = self.nExSize - (self.lnExEndPs[-1] - self.lnExStartPs[-1])
		else:
			nLastEx = self.nExSize - (self.lnExEndPs[0] - self.lnExStartPs[0])
		self.bNMD = (nLastEx - nAfterStop > 50)

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
		return self.sORFSeq
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

	def getNMD(self):
		return self.bNMD
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
		#End for
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
		#End for
		cRefSeq.putCDSSeq(sCDSSeq)
		lProcRefSeq.append(cRefSeq)
	#End for
	return lProcRefSeq

def takecareNegStrand(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.putTrueEx()
		cRefSeq.putTrueCDS()
		lProcCRefSeq.append(cRefSeq)
	#End for
	return lProcCRefSeq

def genORFSeqs(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.parseORF()
		lProcCRefSeq.append(cRefSeq)
	#End for
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
	#End for
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
	#End for
	return listFiltered

def fillNMD(listCRefSeq):
	lProcCRefSeq = []
	for cRefSeq in listCRefSeq:
		cRefSeq.checkNMD()
		lProcCRefSeq.append(cRefSeq)
	#End for
	return lProcCRefSeq

#Mission2 Functions
def parseWholeRefFlat(sFileName):
	wholeCRefSeq = []
	fileH = open(sFileName, 'r')
	for sLine in fileH:
		cRefSeq = RefSeq()
		cRefSeq.parse_refflat_line(sLine)
		wholeCRefSeq.append(cRefSeq)
	#End for
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
	#End for
	#outH = open("Excluded", 'w')
	#print("not NM_", file=outH)
	#[print(sRefID, file=outH) for sRefID in excludedRefID]
	#print("\n\nnot Chrom Num", file=outH)
	#[print(sChromID, file=outH) for sChromID in excludedChromID]
	#outH.close()
	return filteredList

def parseUniqueNM(listCRefSeq):
	refDict = defaultdict(list)
	for cRefSeq in listCRefSeq:
		sRefID = cRefSeq.getRefID()
		refDict[sRefID].append(cRefSeq)
	#End for
	uniqueList = [cRefSeq for cRefSeq in listCRefSeq if len(refDict[cRefSeq.getRefID()]) == 1]
	overlapList = set([cRefSeq.getRefID() for cRefSeq in listCRefSeq if len(refDict[cRefSeq.getRefID()]) != 1])

	#outH = open("NotUnique", 'w')
	#print("not Unique", file=outH)
	#[print(sRefID, file=outH) for sRefID in overlapList]
	#outH.close()

	return uniqueList

def sortByRefSeqID(listCRefSeq):
	listCRefSeq.sort(key = (lambda cRefSeq: int(cRefSeq.getRefID().split("_")[-1])))
	return listCRefSeq

def main_3():
	sFileName = "../data/refFlat.txt"
	#sOutFile = "Mission3_HJeon.txt"
	#hOutF = open(sOutFile, 'w')
	listCRefSeq = parseWholeRefFlat(sFileName)
	#print("Answer 1: {}".format(len(listCRefSeq)), file=hOutF)
	listCRefSeq = parseNM(listCRefSeq)
	#print("Answer 2: {}".format(len(listCRefSeq)), file=hOutF)
	listCRefSeq = parseUniqueNM(listCRefSeq)
	#print("Answer 3: {}".format(len(listCRefSeq)), file=hOutF)

	chromDict = {'hg38': ([str(i) for i in range(1,23)]+['X','Y'])}
	listCRefSeq = genExSeqsPool(chromDict, listCRefSeq)
	listCRefSeq = takecareNegStrand(listCRefSeq)
	listCRefSeq = genORFSeqs(listCRefSeq)

	listCRefSeq = onlyValidORF(listCRefSeq)
	#print("Answer 4: {}".format(len(listCRefSeq)), file=hOutF)

	listCRefSeq = sortByRefSeqID(listCRefSeq)

	listCRefSeq = isoformFiltering(listCRefSeq)
	#print("Answer 5: {}".format(len(listCRefSeq)), file=hOutF)

	#hOutF.close()

	listCRefSeq = sortByRefSeqID(listCRefSeq) #Just redundantly sorting...

	return listCRefSeq

	#listRefID = [cRefSeq.getRefID() for cRefSeq in listCRefSeq]
	#listGeneSymbol = [cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]
	#list5UTRSize = [cRefSeq.get5UTRSize() for cRefSeq in listCRefSeq]
	#listORFSize = [cRefSeq.getORFSize() for cRefSeq in listCRefSeq]
	#list3UTRSize = [cRefSeq.get3UTRSize() for cRefSeq in listCRefSeq]
"""
	dfRefSeq = DataFrame({"Ref Seq ID": listRefID, "Gene Symbol": listGeneSymbol,
		"5'UTR Size": list5UTRSize, "ORF Size": listORFSize, "3'UTR Size": list3UTRSize})
	dfRefSeq.to_excel("Mission3_HJeon.xlsx", sheet_name="sheet1", index=False)
	dfRefSeq.to_csv("Mission3_HJeon.csv", sep="\t")

	#Double checking
	for cRefSeq in listCRefSeq:
		sORFSeq = cRefSeq.getORFSeq()
		assert("N" not in sORFSeq)
		assert(cRefSeq.getCDSSeq() == cRefSeq.getORFSeq())
		assert(sORFSeq[:3] == "ATG")
		assert(sORFSeq[-3:] in ["TAA", "TGA", "TAG"])
		assert(len(sORFSeq) % 3 == 0)
		lCodons = []
		for i in range(3, len(sORFSeq) - 3, 3):
			lCodons.append(sORFSeq[i:i+3])
		for stopC in ["TAA", "TGA", "TAG"]:
			assert(stopC not in lCodons)
		assert(cRefSeq.getORFSize() + cRefSeq.get5UTRSize() + cRefSeq.get3UTRSize() == cRefSeq.getExSize())
		print("{} have valid ORF".format(cRefSeq.getRefID()))
	#End for
"""
	#Check NMD
"""
	NMDlist = fillNMD(listCRefSeq)
	NMDlist = [cRefSeq for cRefSeq in NMDlist if cRefSeq.getNMD()]
	NMDfile = "NMD_List.txt"
	hNMDF = open(NMDfile, 'w')
	print(len(NMDlist), file=hNMDF)
	print("\n", file=hNMDF)
	for cRefSeq in NMDlist:
		print(cRefSeq.getRefID(), file=hNMDF)
		#print(cRefSeq.getExSeq(), file=hNMDF)
		#print(cRefSeq.getORFSeq(), file=hNMDF)
		if cRefSeq.getStrand() == "+":
			nlast2ExLen = cRefSeq.getExEndPs()[-1] - cRefSeq.getExStartPs()[-1] + cRefSeq.getExEndPs()[-2] - cRefSeq.getExStartPs()[-2]
		else:
			nlast2ExLen = cRefSeq.getExEndPs()[0] - cRefSeq.getExStartPs()[0] + cRefSeq.getExEndPs()[1] - cRefSeq.getExStartPs()[1]
		print(cRefSeq.getExSize() - nlast2ExLen - (cRefSeq.getExSize() - cRefSeq.get3UTRSize()), file=hNMDF)
		print("\n", file=hNMDF)
	hNMDF.close()
"""
	#Checking Seqs
"""
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

WHOLE_SIZE = 0;
DOWN_SIZE = 0;
NOT_DOWN_SIZE = 0;

class Motif:
	def __init__(self):
		self.sMotif = ""
		self.fPValue = 0.
		self.fRelRisk = 0.
		self.mtxContingency = np.zeros((2, 2))
		self.nSizeDownWMotif = 0
		self.nSizeDownWNMotif = 0
		self.nSizeNDownWMotif = 0
		self.nSizeNDownWNMotif = 0
		self.dict_miRNA_Info = {}

	#Put fuctions
	def putMotif(self, sMotif):
		self.sMotif = sMotif
	def putPValue(self, fPValue):
		self.fPValue = fPValue
	def putRelRisk(self, fRelRisk):
		self.fRelRisk = fRelRisk
	def putContingency(self, mtxContingency):
		self.mtxContingency = mtxContingency
		self.nSizeDownWMotif = mtxContingency[0, 0]
		self.nSizeDownWNMotif = mtxContingency[0, 1]
		self.nSizeNDownWMotif = mtxContingency[1, 0]
		self.nSizeNDownWNMotif = mtxContingency[1, 1]
	def addContingency(self, i, j):
		self.mtxContingency[i, j] += 1
		self.nSizeDownWMotif = self.mtxContingency[0, 0]
		self.nSizeDownWNMotif = self.mtxContingency[0, 1]
		self.nSizeNDownWMotif = self.mtxContingency[1, 0]
		self.nSizeNDownWNMotif = self.mtxContingency[1, 1]
	def fillNotMotifContingency(self, nSizeDown, nSizeNotDown):
		self.mtxContingency[0, 1] = nSizeDown - self.mtxContingency[0, 0]
		self.mtxContingency[1, 1] = nSizeNotDown - self.mtxContingency[1, 0]
		self.nSizeDownWMotif = self.mtxContingency[0, 0]
		self.nSizeDownWNMotif = self.mtxContingency[0, 1]
		self.nSizeNDownWMotif = self.mtxContingency[1, 0]
		self.nSizeNDownWNMotif = self.mtxContingency[1, 1]
	def put_miRNA_Info(self, sName, sType):
		self.dict_miRNA_Info[sName] = sType

	def getParsedTypeInfo(self):
		sTypeInfo = ""

		for sMiRNA, sType in self.dict_miRNA_Info.items():
			sTypeInfo += sType + ": " + sMiRNA + ","

		return sTypeInfo.strip(',')

	#Get fucntions
	def getMotif(self):
		return self.sMotif
	def getPValue(self):
		return self.fPValue
	def getRelRisk(self):
		return self.fRelRisk
	def getSizeDownWMotif(self):
		return self.nSizeDownWMotif
	def getSizeDownWNMotif(self):
		return self.nSizeDownWNMotif
	def getSizeNDownWMotif(self):
		return self.nSizeNDownWMotif
	def getSizeNDownWNMotif(self):
		return self.nSizeNDownWNMotif
	def getContingency(self):
		return self.mtxContingency
	def get_miRNA_Info(self):
		return self.dict_miRNA_Info

def parseRegData(sRegDataFile):
	dictRegData = {}
	hRegDataF = open(sRegDataFile, 'r')
	sRegData = hRegDataF.readlines()
	hRegDataF.close()

	for line in sRegData:
		flds = line.split('\t')
		assert(len(flds)==2), "{}: Does not have proper format".format(line)
		sGeneSymbol = flds[0].strip().upper()
		fReg = float(flds[1].strip())
		dictRegData[sGeneSymbol] = fReg
	#End for
	return dictRegData

def genEveryMotif(nMotifLen, listCRefSeq, sRegion):
	setEveryMotif = set([])
	if sRegion == "3UTR":
		for cRefSeq in listCRefSeq:
			s3UTRSeq = cRefSeq.get3UTRSeq()
			for i in range(len(s3UTRSeq) - nMotifLen + 1):
				sMotif = s3UTRSeq[i : i+nMotifLen]
				setEveryMotif.add(sMotif)

	if sRegion == "ORF":
		for cRefSeq in listCRefSeq:
			sORFSeq = cRefSeq.getORFSeq()
			for i in range(len(sORFSeq) - nMotifLen + 1):
				sMotif = sORFSeq[i : i+nMotifLen]
				setEveryMotif.add(sMotif)

	if sRegion == "5UTR":
		for cRefSeq in listCRefSeq:
			s5UTRSeq = cRefSeq.get5UTRSeq()
			for i in range(len(s5UTRSeq) - nMotifLen + 1):
				sMotif = s5UTRSeq[i : i+nMotifLen]
				setEveryMotif.add(sMotif)

	return list(setEveryMotif)

def initMotifClasses(listEveryMotif):
	dictCMotif = {}
	for sMotif in listEveryMotif:
		cMotif = Motif()
		cMotif.putMotif(sMotif)
		dictCMotif[sMotif] = cMotif

	return dictCMotif

def fillMotifClasses(dictRegData, dictCMotif, listCRefSeq, sRegion):
	#sNotInReg = "Not_in_reg_opt.txt"
	#hNotInReg = open(sNotInReg, 'w')

	nMotifLen = len(list(dictCMotif.keys())[0])
	print("Motif Length: {}".format(nMotifLen))
	
	listCRefSeq = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getGeneSymbol().upper() in dictRegData.keys()]
	listDownCRef = [cRefSeq for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] < -0.5]
	listNotDownCRef = [cRefSeq for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] >= -0.5]
	assert(len(listDownCRef) + len(listNotDownCRef) == len(listCRefSeq))

	#For checking
	global WHOLE_SIZE
	global DOWN_SIZE
	global NOT_DOWN_SIZE

	WHOLE_SIZE = len(listCRefSeq)
	DOWN_SIZE = len(listDownCRef)
	NOT_DOWN_SIZE = len(listNotDownCRef)

	nDownSize = len(listDownCRef)
	nNotDownSize = len(listNotDownCRef)

	print("Down Size: {}".format(nDownSize))
	print("Not Down Size: {}".format(nNotDownSize))

	if sRegion == "3UTR":
		for cRefSeq in listDownCRef:
			s3UTR = cRefSeq.get3UTRSeq()
			n3UTRSize = len(s3UTR)
			setSeenMotif = set([])

			for i in range(n3UTRSize - nMotifLen + 1):
				sMotif = s3UTR[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(0, 0)
					dictCMotif[sMotif] = cMotif

		for cRefSeq in listNotDownCRef:
			s3UTR = cRefSeq.get3UTRSeq()
			n3UTRSize = len(s3UTR)
			setSeenMotif = set([])

			for i in range(n3UTRSize - nMotifLen + 1):
				sMotif = s3UTR[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(1, 0)
					dictCMotif[sMotif] = cMotif

		for key, val in dictCMotif.items():
			cMotif = val
			cMotif.fillNotMotifContingency(nDownSize, nNotDownSize)
			dictCMotif[key] = cMotif

	if sRegion == "ORF":
		for cRefSeq in listDownCRef:
			sORF = cRefSeq.getORFSeq()
			nORFSize = len(sORF)
			setSeenMotif = set([])

			for i in range(nORFSize - nMotifLen + 1):
				sMotif = sORF[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(0, 0)
					dictCMotif[sMotif] = cMotif

		for cRefSeq in listNotDownCRef:
			sORF = cRefSeq.getORFSeq()
			nORFSize = len(sORF)
			setSeenMotif = set([])

			for i in range(nORFSize - nMotifLen + 1):
				sMotif = sORF[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(1, 0)
					dictCMotif[sMotif] = cMotif

		for key, val in dictCMotif.items():
			cMotif = val
			cMotif.fillNotMotifContingency(nDownSize, nNotDownSize)
			dictCMotif[key] = cMotif

	if sRegion == "5UTR":
		for cRefSeq in listDownCRef:
			s5UTR = cRefSeq.get5UTRSeq()
			n5UTRSize = len(s5UTR)
			setSeenMotif = set([])

			for i in range(n5UTRSize - nMotifLen + 1):
				sMotif = s5UTR[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(0, 0)
					dictCMotif[sMotif] = cMotif

		for cRefSeq in listNotDownCRef:
			s5UTR = cRefSeq.get5UTRSeq()
			n5UTRSize = len(s5UTR)
			setSeenMotif = set([])

			for i in range(n5UTRSize - nMotifLen + 1):
				sMotif = s5UTR[i : i+nMotifLen]
				if sMotif not in setSeenMotif:
					setSeenMotif.add(sMotif)
					cMotif = dictCMotif[sMotif]
					cMotif.addContingency(1, 0)
					dictCMotif[sMotif] = cMotif

		for key, val in dictCMotif.items():
			cMotif = val
			cMotif.fillNotMotifContingency(nDownSize, nNotDownSize)
			dictCMotif[key] = cMotif

	#hNotInReg.close()
	return list(dictCMotif.values())

def fillStatistics(listCMotif):
	listStatMotif = []
	for cMotif in listCMotif:
		mtxContingency = cMotif.getContingency()
		OddsRatio, PValue = stats.fisher_exact(mtxContingency)
		fRelRisk = calculateRelRisk(mtxContingency)
		cMotif.putPValue(PValue)
		cMotif.putRelRisk(fRelRisk)
		listStatMotif.append(cMotif)

	return listStatMotif

def calculateRelRisk(mtxContingency):
	fN1 = mtxContingency[0, 0]
	fN2 = mtxContingency[0, 1]
	fN3 = mtxContingency[1, 0]
	fN4 = mtxContingency[1, 1]

	fA = fN1 / (fN1 + fN2)
	fB = fN3 / (fN3 + fN4)

	fRelRisk = fA / fB

	return fRelRisk

def significantCMotifs(listCMotif):
	sigCMotif = [cMotif for cMotif in listCMotif if cMotif.getRelRisk() > 1.0]
	sigCMotif.sort(key = (lambda cMotif: cMotif.getPValue()))

	return sigCMotif



def main_4(sRegDataFile, listCRefSeq, sRegion):
	nMotifLen = 7

	#listCRefSeq = main_3()
	assert(len(listCRefSeq) == 19076)
	dictRegData = parseRegData(sRegDataFile)

	#sInRegNotInRef = "In_reg_not_in_ref_faster.txt"
	#hInRegNotInRef = open(sInRegNotInRef, 'w')
	#listInRegNotInRef = list(set(list(dictRegData.keys())) - set([cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]))
	#print("\n".join(listInRegNotInRef), file=hInRegNotInRef)
	#hInRegNotInRef.close()

	listEveryMotif = genEveryMotif(nMotifLen, listCRefSeq, sRegion)
	print("Number of Every Motifs: {}".format(len(listEveryMotif)))
	dictCMotif = initMotifClasses(listEveryMotif)
	print("Length of cMotif List: {}".format(len(list(dictCMotif.keys()))))
	assert(len(listEveryMotif) == len(dictCMotif.keys()))
	listCMotif = fillMotifClasses(dictRegData, dictCMotif, listCRefSeq, sRegion)
	listCMotif = fillStatistics(listCMotif)

	#Filter RelRisk > 1.0
	#Sort with P Value
	listCMotif = significantCMotifs(listCMotif)

	return listCMotif

	#Write to Excel
	"""
	listMotifSeq = [cMotif.getMotif() for cMotif in listCMotif]
	listPValue = [cMotif.getPValue() for cMotif in listCMotif]
	listMotif_Down = [cMotif.getSizeDownWMotif() for cMotif in listCMotif]
	listNotMotif_Down = [cMotif.getSizeDownWNMotif() for cMotif in listCMotif]
	listMotif_NotDown = [cMotif.getSizeNDownWMotif() for cMotif in listCMotif]
	listNotMotif_NotDown = [cMotif.getSizeNDownWNMotif() for cMotif in listCMotif]
	listRelRisk = [cMotif.getRelRisk() for cMotif in listCMotif]

	print("Down Size: {}".format(DOWN_SIZE))
	arrayDown = np.array(listMotif_Down) + np.array(listNotMotif_Down)
	assert(np.all(arrayDown == DOWN_SIZE))

	print("Not Down Size: {}".format(NOT_DOWN_SIZE))
	arrayNotDown = np.array(listMotif_NotDown) + np.array(listNotMotif_NotDown)
	assert(np.all(arrayNotDown == NOT_DOWN_SIZE))

	print("Whole Size: {}".format(WHOLE_SIZE))
	arrayWhole = arrayDown + arrayNotDown
	assert(np.all(arrayWhole == WHOLE_SIZE))

	dfMotif = DataFrame({"Motif": listMotifSeq, "P_Value": listPValue,
		"Motif_Down": listMotif_Down, "NotMotif_Down": listNotMotif_Down,
		"Motif_NotDown": listMotif_NotDown, "NotMotif_NotDown": listNotMotif_NotDown,
		"Relative_Risk": listRelRisk})

	dfMotif.to_excel("Mission4_HJeon.xlsx", sheet_name="sheet1", index=False)
	"""	

def gen_dict_miRNA(listCMiRNA):
	dictCMiRNA = {}

	for cMiRNA in listCMiRNA:
		sName = cMiRNA.getName()
		dictCMiRNA[sName] = cMiRNA

	return dictCMiRNA

def check_miRNA_extension(cMiRNA, listCRefSeq):
	listExtLen = []
	sRNASeq = cMiRNA.getSeq()
	for cRefSeq in listCRefSeq:
		sExSeq = cRefSeq.getORFSeq()
		cnt = 7
		sMotif = sRNASeq[-cnt:]
		if sMotif in sExSeq:
			while sMotif in sExSeq:
				cnt += 1
				sMotif = sRNASeq[-cnt:]
			listExtLen.append(cnt)
	return listExtLen

def check_target_repeats(cMiRNA, listCRefSeq):
	listRepeatLen = []
	sRNASeq = cMiRNA.getSeq()[1:8]
	for cRefSeq in listCRefSeq:
		sORFSeq = cRefSeq.getORFSeq()
		nRepeat_Count = 0
		if sRNASeq in sORFSeq:
			for i in range(len(sORFSeq)):
				sSub = sORFSeq[i:i+7]
				if sRNASeq == sSub:
					nRepeat_Count += 1
			listRepeatLen.append(nRepeat_Count)

	return listRepeatLen 

def check_correlation(cMiRNA, listCRefSeq):
	corr_cnt = 0
	sRNASeq = cMiRNA.getSeq()[-7:]
	for cRefSeq in listCRefSeq:
		if sRNASeq in cRefSeq.get3UTRSeq() and sRNASeq in cRefSeq.getORFSeq():
			corr_cnt += 1

	return corr_cnt

def main_5():
	s_miRNA_Name = "miR-9-5p"
	sRegDataFile = "../data/Mission5_Dataset1.txt"

	pickle_off = open("../data/Mission3.pickle", "rb")
	listCRefSeq = pickle.load(pickle_off)
	pickle_off.close()
	sMiRNAFile = "../data/mature.fa"
	listCMiRNA = fill_miRNA(sMiRNAFile)
	dictCMiRNA = gen_dict_miRNA(listCMiRNA)
	dictRegData = parseRegData(sRegDataFile)
	dictRegData = {k: v for k, v in dictRegData.items() if v < -0.5}
	listCRefSeq = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getGeneSymbol().upper() in dictRegData]
	#listCMotif = main_4(sRegDataFile, listCRefSeq, sRegion)
	#cMotif = [cMotif for cMotif in listCMotif if cMotif.getMotif() == sMotif][0]
	cMiRNA = dictCMiRNA[s_miRNA_Name]
	nNumCorr = check_correlation(cMiRNA, listCRefSeq)
	print(nNumCorr)

if __name__ == "__main__":
	rtime = time()
	main_5()
	rtime = time() - rtime
	print("{} seconds".format(rtime))