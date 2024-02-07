#!/usr/bin/env python
import glob
import os
from time import sleep

def IntersectBed(sPathFile):
    sOutFile=sPathFile.split("/")[-1]
    os.system("bedtools intersect -wa -a /ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/Sorted_Interval_paz.bed -b "+sPathFile+" -loj > /ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/Intersected_paz_Cancergroup/Intersected_"+sOutFile)
   # sleep(1)
    os.system("bedtools intersect -wa -c -a /ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/Sorted_Interval_paz.bed -b "+sPathFile+" -loj > /ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/IntersectedCount_paz_Cancergroup/IntersectedCount_"+sOutFile)
    #sleep(1)






if __name__=="__main__":
#    lFilelists=glob.glob("/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/BarcodeGroup/*.bed")
    lFilelists=glob.glob("/ahg/regevdata/projects/ICA_Lung/Wooseung/ICGC/PCAWG/DIG/CancerGroup/Bed/*.bed")
 
#lFilelists=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/BarcodeGroup/LAML-KR.bed"]
    for sPathFile in lFilelists:
        IntersectBed(sPathFile)







