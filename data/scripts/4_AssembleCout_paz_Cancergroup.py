#!/usr/bin/env python
import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str, default=None)
parser.add_argument('--subsampled', type=bool, default=False)

config = parser.parse_args()
cancer_types = config.cancer_types
subsampled = config.subsampled

current_dir = os.path.dirname(os.path.abspath(__file__))

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

def helper(cancer_type, current_dir, dirname, subsampled):
    if subsampled:
        sOutname=f"{cancer_type}_{os.path.basename(dirname)}.txt"
    else:
        sOutname=f"{cancer_type}.txt"

    lFilelist=glob.glob(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{dirname}/IntersectedCount_paz_Cancergroup/*.bed")
    #lFilelist=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/IntersectedCount/IntersectedCount_LAML-KR.bed"]
    dTumorCount=dict()
    for sPathFile in lFilelist:
        sFile=sPathFile.split("/")[-1]
        sTumorname=sFile.split(".")[0]
        sTumorBarcode=sTumorname.split("_")[1]
        dTumorCount[sTumorBarcode]=ParsetoDict(sPathFile)

    fp=open(f"{current_dir}/../Sorted_Interval_paz.bed")
    lKeylist=[]
    for sLine in fp.readlines():
        sLine=sLine.strip()
        t=sLine.split("\t")
        (sChr,sStart,sEnd,sSeqname)=(t[0],t[1],t[2],t[3])
        sKey=sChr+"_"+sStart+"_"+sEnd+"_"+sSeqname
        lKeylist.append(sKey)
    fout=open(f"{current_dir}/../processed_data/{sOutname}", "w")
    fout.write("Chr\tStart\tEnd\tSeqname\t")
    lTumorKey=dTumorCount.keys()
    fout.write("{0}\n".format("\t".join(lTumorKey)))
    for sKey in lKeylist:
        lTemporlist=[]
        fout.write("{0}\t".format("\t".join(sKey.split("_"))))
        for sTumor in lTumorKey:
            lTemporlist.append(dTumorCount[sTumor][sKey])
        fout.write("{0}\n".format("\t".join(map(str,lTemporlist))))

if __name__=="__main__":
    for cancer_type in cancer_types:
        if subsampled:
            for dirname in os.listdir(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled"):
                helper(cancer_type, current_dir, dirname, subsampled)
        else:
            helper(cancer_type, current_dir, f"{current_dir}/../mutation_data/bed_files/{cancer_type}", subsampled)
        # for dirname in os.listdir(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled"):
        #     # tStart=time.ctime()
        #     sOutname=f"{cancer_type}_{os.path.basename(dirname)}.txt"
        #
        #     lFilelist=glob.glob(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{dirname}/IntersectedCount_paz_Cancergroup/*.bed")
        #     #lFilelist=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/IntersectedCount/IntersectedCount_LAML-KR.bed"]
        #     dTumorCount=dict()
        #     for sPathFile in lFilelist:
        #         sFile=sPathFile.split("/")[-1]
        #         sTumorname=sFile.split(".")[0]
        #         sTumorBarcode=sTumorname.split("_")[1]
        #         dTumorCount[sTumorBarcode]=ParsetoDict(sPathFile)
        #
        #     fp=open(f"{current_dir}/../Sorted_Interval_paz.bed")
        #     lKeylist=[]
        #     for sLine in fp.readlines():
        #         sLine=sLine.strip()
        #         t=sLine.split("\t")
        #         (sChr,sStart,sEnd,sSeqname)=(t[0],t[1],t[2],t[3])
        #         sKey=sChr+"_"+sStart+"_"+sEnd+"_"+sSeqname
        #         lKeylist.append(sKey)
        #     fout=open(f"{current_dir}/../processed_data/{sOutname}", "w")
        #     fout.write("Chr\tStart\tEnd\tSeqname\t")
        #     lTumorKey=dTumorCount.keys()
        #     fout.write("{0}\n".format("\t".join(lTumorKey)))
        #     for sKey in lKeylist:
        #         lTemporlist=[]
        #         fout.write("{0}\t".format("\t".join(sKey.split("_"))))
        #         for sTumor in lTumorKey:
        #             lTemporlist.append(dTumorCount[sTumor][sKey])
        #         fout.write("{0}\n".format("\t".join(map(str,lTemporlist))))

    # print(tStart)
    # print(time.ctime())







