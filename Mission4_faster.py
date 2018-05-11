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
	dictCMotif = {}
	for sMotif in listEveryMotif:
		cMotif = Motif()
		cMotif.putMotif(sMotif)
		dictCMotif[sMotif] = cMotif

	return dictCMotif

def fillMotifClasses(dictRegData, dictCMotif, listCRefSeq):
	#sNotInReg = "Not_in_reg_opt.txt"
	#hNotInReg = open(sNotInReg, 'w')

	nMotifLen = len(list(dictCMotif.keys())[0])
	print("Motif Length: {}".format(nMotifLen))
	
	listCRefSeq = [cRefSeq for cRefSeq in listCRefSeq if cRefSeq.getGeneSymbol().upper() in dictRegData.keys()]
	listDownCRef = [cRefSeq for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] < -0.5]
	listNotDownCRef = [cRefSeq for cRefSeq in listCRefSeq if dictRegData[cRefSeq.getGeneSymbol().upper()] >= -0.5]

	nDownSize = len(listDownCRef)
	nNotDownSize = len(listNotDownCRef)

	print("Down Size: {}".format(nDownSize))
	print("Not Down Size: {}".format(nNotDownSize))

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

def main():
	sRegDataFile = "../data/Mission4_Dataset.txt"
	nMotifLen = 7

	listCRefSeq = Mission3_Lib.main()
	assert(len(listCRefSeq) == 19076)
	dictRegData = parseRegData(sRegDataFile)

	sInRegNotInRef = "In_reg_not_in_ref_faster.txt"
	hInRegNotInRef = open(sInRegNotInRef, 'w')
	listInRegNotInRef = list(set(list(dictRegData.keys())) - set([cRefSeq.getGeneSymbol() for cRefSeq in listCRefSeq]))
	print("\n".join(listInRegNotInRef), file=hInRegNotInRef)
	hInRegNotInRef.close()

	listEveryMotif = genEveryMotif(nMotifLen, listCRefSeq)
	print("Number of Every Motifs: {}".format(len(listEveryMotif)))
	dictCMotif = initMotifClasses(listEveryMotif)
	print("Length of cMotif List: {}".format(len(list(dictCMotif.keys()))))
	assert(len(listEveryMotif) == len(dictCMotif.keys()))
	listCMotif = fillMotifClasses(dictRegData, dictCMotif, listCRefSeq)
	listCMotif = fillStatistics(listCMotif)

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

	dfMotif = DataFrame({"Motif": listMotifSeq, "P_Value": listPValue,
		"Motif_Down": listMotif_Down, "NotMotif_Down": listNotMotif_Down,
		"Motif_NotDown": listMotif_NotDown, "NotMotif_NotDown": listNotMotif_NotDown,
		"Relative_Risk": listRelRisk})

	dfMotif.to_excel("Mission4_faster.xlsx", sheet_name="sheet1", index=False)

if __name__ == "__main__":
	rtime = time()
	main()
	rtime = time() - rtime
	print("{} seconds".format(rtime))