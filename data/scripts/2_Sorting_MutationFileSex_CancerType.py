#!/usr/bin/env python

import glob 
import time
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str, default=None)
config = parser.parse_args()
cancer_types = config.cancer_types

tStartTime = (time.ctime())

current_dir = os.path.dirname(os.path.abspath(__file__))
class cVariant:
	def __init__(self): 
		self.sChr = 'NULL'
		self.nChr = 0
		self.nStart = 0
		self.nEnd = 0
		self.sStrand = 'NULL'
		self.donor_id = 'NULL'

	def parse_line(self, sReadLine):
		sList = sReadLine.split("\t")
		sChr = sList[1].replace("chr","")
		self.sChr = sChr
		self.nStart = int(sList[2])
		self.nEnd = int(sList[3])
		self.sGenename = sList[0]
		self.donor_id = sList[42]
		nChr = sChr.replace("chr","")
		if nChr == "X":
			nChr = 23
		elif nChr == "Y":
			nChr = 24
		else:
			nChr = int(nChr)
		self.nChr = nChr

def SortingMutation(sPathFile):
	fp=open(sPathFile)
	sFile=sPathFile.split("/")[-1]
	sOutname=sFile.split("_")[0]
	lMutationlist=[]
	#fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		cMutationcVariant=cVariant()
		cMutationcVariant.parse_line(sLine)
		lMutationlist.append(cMutationcVariant)
	lMutationlist.sort(key=lambda x:x.nStart)
	lMutationlist.sort(key=lambda x:x.nChr)
	fout = open(f"{current_dir}/../mutation_data/bed_files/{sOutname}.bed", "w")
	for cMutationcVariant in lMutationlist:
		fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(cMutationcVariant.sChr,
													  cMutationcVariant.nStart,
													  cMutationcVariant.nEnd,
													  cMutationcVariant.sGenename,
													  cMutationcVariant.donor_id))

if __name__ == "__main__":
	tStartTime=(time.ctime())
	lFilelist = []
	for cancer_type in cancer_types:
		pattern = current_dir + f"/../mutation_data/{cancer_type}_SNV_with_SEX.txt"
		lFilelist = lFilelist + glob.glob(pattern)
	for sFile in lFilelist:
		SortingMutation(sFile)
	print("Start Time")
	print(tStartTime)
	print("End Time")
	print(time.ctime())























			








