#!/usr/bin/env python

import glob 
import re 
import time 
tStartTime= (time.ctime()) 


class cVariant: 
 
	def __init__(self): 
		self.sChr =     'NULL' 
		self.nChr=	0
		self.nStart =  0
		self.nEnd =    0 
		self.sStrand=        'NULL' 

	def parse_line(self,sReadLine):

		sList=sReadLine.split("\t")
		sChr=sList[1].replace("chr","") 
		self.sChr=      sChr
		self.nStart=     int(sList[2])
		self.nEnd=     int(sList[3])
		self.sGenename=      sList[0]      

		
		nChr=sChr.replace("chr","")
		if nChr=="X":
			nChr=23
		elif nChr=="Y":
			nChr=24
		else:
			nChr=int(nChr)
		self.nChr=nChr


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
	fout=open("/ahg/regevdata/projects/ICA_Lung/Wooseung/ICGC/PCAWG/DIG/CancerGroup/Bed/"+sOutname+".bed","w")
	for cMutationcVariant in lMutationlist:
		fout.write("{0}\t{1}\t{2}\t{3}\n".format(cMutationcVariant.sChr,cMutationcVariant.nStart,cMutationcVariant.nEnd,cMutationcVariant.sGenename))










if __name__=="__main__":

	tStartTime= (time.ctime()) 


	lFilelist=glob.glob("/ahg/regevdata/projects/ICA_Lung/Wooseung/ICGC/PCAWG/DIG/CancerGroup/*SNV_with_SEX.txt")
	#lFilelist=["/ahg/regevdata/projects/ICA_Lung/Wooseung/ICGC/PCAWG/DIG/BarcodeGroup/SNV_without_SEX/LAML-KR_SNV_without_SEX.txt"]

	for sFile in lFilelist:
		SortingMutation(sFile)


	print("Start Time")
	print(tStartTime)
	print("End Time")
	print(time.ctime())























			








