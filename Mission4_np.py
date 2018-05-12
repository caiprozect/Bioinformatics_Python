from time import time
from scipy import stats
from pandas import DataFrame
import numpy as np
import pickle

WHOLE_SIZE = 0;
DOWN_SIZE = 0;
NOT_DOWN_SIZE = 0;

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

def genEveryMotif(nMotifLen, listCRefSeq):
	setEveryMotif = set([])
	for cRefSeq in listCRefSeq:
		s3UTRSeq = cRefSeq.get3UTRSeq()
		for i in range(len(s3UTRSeq) - nMotifLen + 1):
			sMotif = s3UTRSeq[i : i+nMotifLen]
			setEveryMotif.add(sMotif)

	return list(setEveryMotif)

def initMotifClasses(listEveryMotif):
	listCMotif = []
	for sMotif in listEveryMotif:
		cMotif = Motif()
		cMotif.putMotif(sMotif)
		listCMotif.append(cMotif)

	return listCMotif

def fillMotifClasses(dictRegData, listCMotif, listCRefSeq):
	#sNotInReg = "Not_in_reg_opt.txt"
	#hNotInReg = open(sNotInReg, 'w')

	listFilledMotif = []

	nMotifLen = len(listCMotif[0].getMotif())
	print("Motif Length: {}".format(nMotifLen))
	
	listCRefSeq = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getGeneSymbol().upper() in dictRegData.keys()]
	listDownCRef = [cRefSeq.get3UTRSeq() for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] < -0.5]
	listNotDownCRef = [cRefSeq.get3UTRSeq() for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] >= -0.5]
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

	for cMotif in listCMotif:
		mtxContingency = np.zeros((2,2))
		sMotif = cMotif.getMotif()
		mtxContingency[0,0] = np.sum(np.char.find(listDownCRef, sMotif) != -1)
		mtxContingency[0,1] = np.sum(np.char.find(listDownCRef, sMotif) == -1)
		mtxContingency[1,0] = np.sum(np.char.find(listNotDownCRef, sMotif) != -1)
		mtxContingency[1,1] = np.sum(np.char.find(listDownCRef, sMotif) == -1)
		cMotif.putContingency(mtxContingency)
		OddsRatio, PValue = stats.fisher_exact(mtxContingency)
		fRelRisk = calculateRelRisk(mtxContingency)
		cMotif.putPValue(PValue)
		cMotif.putRelRisk(fRelRisk)
		listFilledMotif.append(cMotif)

	#hNotInReg.close()
	return listFilledMotif

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

def main():
	sRegDataFile = "../data/Mission4_Dataset.txt"
	nMotifLen = 7

	pickle_off = open("../data/Mission3.pickle", "rb")
	listCRefSeq = pickle.load(pickle_off)
	pickle_off.close()
	assert(len(listCRefSeq) == 19076)
	dictRegData = parseRegData(sRegDataFile)

	#sInRegNotInRef = "In_reg_not_in_ref_np.txt"
	#hInRegNotInRef = open(sInRegNotInRef, 'w')
	#listInRegNotInRef = list(set(list(dictRegData.keys())) - set([cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]))
	#print("\n".join(listInRegNotInRef), file=hInRegNotInRef)
	#hInRegNotInRef.close()

	listEveryMotif = genEveryMotif(nMotifLen, listCRefSeq)
	print("Number of Every Motifs: {}".format(len(listEveryMotif)))
	listCMotif = initMotifClasses(listEveryMotif)
	print("Length of cMotif List: {}".format(len(listCMotif)))
	assert(len(listEveryMotif) == len(listCMotif))
	listCMotif = fillMotifClasses(dictRegData, listCMotif, listCRefSeq)

	#Filter RelRisk > 1.0
	#Sort with P Value
	listCMotif = significantCMotifs(listCMotif)

	#Write to Excel
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

	dfMotif.to_excel("Mission4_np.xlsx", sheet_name="sheet1", index=False)

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("{} seconds".format(rtime))