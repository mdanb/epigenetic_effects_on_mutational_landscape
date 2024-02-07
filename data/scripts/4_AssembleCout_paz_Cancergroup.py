#!/usr/bin/env python
import glob
import time

def ParsetoDict(sPathFile):
    dInterDict=dict()
    fp=open(sPathFile)
    for sLine in fp.readlines():
        sLine=sLine.strip()
        t=sLine.split("\t")
        (sChr,sStart,sEnd,sSeqname,sCount)=(t[0],t[1],t[2],t[3],t[4])
        sKey=sChr+"_"+sStart+"_"+sEnd+"_"+sSeqname
        dInterDict[sKey]=sCount
    return dInterDict




if __name__=="__main__":
    tStart=time.ctime()
    sOutname="PCAWG_Cancergroup_MutationCountperMB.txt"

    lFilelist=glob.glob("/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/IntersectedCount_paz_Cancergroup/*.bed")
    #lFilelist=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/IntersectedCount/IntersectedCount_LAML-KR.bed"]
    dTumorCount=dict()
    for sPathFile in lFilelist:
        sFile=sPathFile.split("/")[-1]
        sTumorname=sFile.split(".")[0]
        sTumorBarcode=sTumorname.split("_")[1]
        dTumorCount[sTumorBarcode]=ParsetoDict(sPathFile)

    fp=open("/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/Sorted_Interval_paz.bed")
    lKeylist=[]
    for sLine in fp.readlines():
        sLine=sLine.strip()
        t=sLine.split("\t")
        (sChr,sStart,sEnd,sSeqname)=(t[0],t[1],t[2],t[3])
        sKey=sChr+"_"+sStart+"_"+sEnd+"_"+sSeqname
        lKeylist.append(sKey)
    fout=open("/ahg/regevdata/projects/ICA_Lung/Wooseung/"+sOutname,"w")
    fout.write("Chr\tStart\tEnd\tSeqname\t")
    lTumorKey=dTumorCount.keys()
    fout.write("{0}\n".format("\t".join(lTumorKey)))
    for sKey in lKeylist:
        lTemporlist=[]
        fout.write("{0}\t".format("\t".join(sKey.split("_"))))
        for sTumor in lTumorKey:
            lTemporlist.append(dTumorCount[sTumor][sKey])
        fout.write("{0}\n".format("\t".join(map(str,lTemporlist))))


    print(tStart)
    print(time.ctime())







