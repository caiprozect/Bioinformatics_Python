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
	listCMotif = []
	for sMotif in listEveryMoitf:
		cMotif = Motif()
		cMotif.putMotif(sMotif)
		listCMotif.append(cMotif)

	return listCMotif

def fillMotifClasses(dictRegData, listCMotif, listCRefSeq):
	listFilledCMotif = []
	for cMotif in listCMotif:
		mtxContingency = np.zeros((2,2))
		sMotif = cMotif.getMotif()
		for cRefSeq in listCRefSeq:
			i = 0
			j = 0
			sGeneSymbol = cRefSeq.getGeneSymbol().upper()
			s3UTRSeq = cRefSeq.get3UTRSeq()
			if sGeneSymbol in dictRegData:
				fReg = dictRegData[sGeneSymbol]
				if (fReg < (-0.05)):
					i = 0
				else:
					i = 1
				checkMotif = s3UTRSeq.find(sMotif)
				if (checkMotif >= 0):
					j = 0
				else:
					j = 1
				mtxContingency[i, j] += 1
			else:
				print("{} Not in Regulation Data".format(sGeneSymbol))
		assert(np.sum(mtxContingency) == len(dictRegData.keys()))
		cMotif.putContingency(mtxContingency)
		fOddsRatio, fPValue = stats.fisher_exact(mtxContingency)
		cMotif.putPValue(fPValue)
		fRelRisk = calculateRelRisk(mtxContingency)
		cMotif.putRelRisk(fRelRisk)
		listFilledCMotif.append(cMotif)
	return listFilledCMotif

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
	listCMotif = initMotifClasses(listEveryMoitf)
	print("Length of cMotif List: {}".format(len(listCMotif)))
	assert(len(listEveryMoitf) == len(listCMotif))
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
	listNotMotif_NotDown = [cMotif.getSizeNMotifWNDown() for cMotif in listCMotif]
	listRelRisk = [cMotif.getRelRisk() for cMotif in listCMotif]

	dfMotif = DataFrame({"Motif": listMotifSeq, "P_Value": listPValue,
		"Motif_Down": listMotif_Down, "NotMotif_Down": listNotMotif_Down,
		"Motif_NotDown": listNotMotif_Down, "NotMotif_NotDown": listNotMotif_NotDown,
		"Relative_Risk": listRelRisk})

	dfMotif.to_excel("Mission4_proto.xlsx", sheet_name="sheet1", index=False)

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("{} seconds".format(rtime))