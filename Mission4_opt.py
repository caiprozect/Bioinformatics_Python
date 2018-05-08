from time import time
from scipy import stats
from pandas import DataFrame
import numpy as np
import Mission3_Lib

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

def initMotifClasses(listEveryMoitf):
	dictCMotif = {}
	for sMotif in listEveryMoitf:
		cMotif = Motif()
		cMotif.putMotif(sMotif)
		dictCMotif[sMotif] = cMotif

	return dictCMotif

def fillMotifClasses(listEveryMotif, dictRegData, dictCMotif, listCRefSeq):
	sNotInReg = "Not_in_reg_opt.txt"
	hNotInReg = open(sNotInReg, 'w')
	nMotifLen = len(list(dictCMotif.keys())[0])
	setEveryMotif = set(listEveryMoitf)
	for cRefSeq in listCRefSeq:
		setSeenMotif = 
		sGeneSymbol = cRefSeq.getGeneSymbol().strip().upper()
		if sGeneSymbol in list(dictRegData.keys()):
			fReg = dictRegData[sGeneSymbol]
			row = 0 if (fReg < -0.5) else 1
			s3UTR = cRefSeq.get3UTRSeq()
			for i in range(len(s3UTR) - nMotifLen + 1):
				sMotif = s3UTR[i:i+nMotifLen]

		else:
			print(sGeneSymbol, file=hNotInReg)
	hNotInReg.close()
	return list(dictCMotif.values())

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

	listCRefSeq = Mission3_Lib.main()
	assert(len(listCRefSeq) == 19076)
	dictRegData = parseRegData(sRegDataFile)
	listEveryMoitf = genEveryMotif(nMotifLen, listCRefSeq)
	print("Number of Every Motifs: {}".format(len(listEveryMoitf)))
	dictCMotif = initMotifClasses(listEveryMoitf)
	print("Length of cMotif List: {}".format(len(listCMotif)))
	assert(len(listEveryMoitf) == len(dictCMotif.keys()))
	listCMotif = fillMotifClasses(listEveryMotif, dictRegData, dictCMotif, listCRefSeq)

	#Filter RelRisk > 1.0
	#Sort with P Value
	listCMotif = significantCMotifs(listCMotif)

	#Write to Excel
	listMotifSeq = [cMotif.getMotif() for cMotif in listCMotif]
	listPValue = [cMotif.getPValue() for cMotif in listCMotif]
	listMotif_Down = [cMotif.getSizeDownWMotif() for cMotif in listCMotif]
	listNotMotif_Down = [cMotif.getSizeDownWNMotif() for cMotif in listCMotif]
	listMotif_NotDown = [cMotif.getSizeNDownWMotif() for cMotif in listCMotif]
	listNotMotif_NotDown = [cMotif.getSizeNMotifWNDown() for cMotif in listCMotif]
	listRelRisk = [cMotif.getRelRisk() for cMotif in listCMotif]

	dfMotif = DataFrame({"Motif": listMotifSeq, "P_Value": listPValue,
		"Motif_Down": listMotif_Down, "NotMotif_Down": listNotMotif_Down,
		"Motif_NotDown": listNotMotif_Down, "NotMotif_NotDown": listNotMotif_NotDown,
		"Relative_Risk": listRelRisk})

	dfMotif.to_excel("Mission4_opt.xlsx", sheet_name="sheet1", index=False)

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("{} seconds".format(rtime))