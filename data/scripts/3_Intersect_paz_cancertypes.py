#!/usr/bin/env python
import glob
import os
import argparse
import pyranges as pr
import numpy as np
import random

random.seed(42)
parser = argparse.ArgumentParser()
parser.add_argument('--cancer_types', nargs="+", type=str, default=None)
parser.add_argument('--decrement_by', type=int, default=None)
config = parser.parse_args()
cancer_types = config.cancer_types
decrement_by = config.decrement_by

current_dir = os.path.dirname(os.path.abspath(__file__))

def IntersectBed(sPathFile, outPath, interval_fp, subsample_idx):
    sOutFile=sPathFile.split("/")[-1]
    sOutFile = sOutFile.split(".")
    sOutFile = f"{sOutFile[0]}_S{subsample_idx}.bed"
    os.system(f"bedtools intersect -wa -a {interval_fp} -b "+sPathFile+ f" -loj > {outPath}/Intersected_paz_Cancergroup/Intersected_"+sOutFile)
   # sleep(1)
    os.system(f"bedtools intersect -wa -c -a {interval_fp} -b "+sPathFile+f" -loj > {outPath}/IntersectedCount_paz_Cancergroup/IntersectedCount_"+sOutFile)
    #sleep(1)

if __name__=="__main__":
    lFilelists = []
    for cancer_type in cancer_types:
        os.makedirs(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/", exist_ok=True)
        lFilelists = lFilelists + glob.glob(f"{current_dir}/{cancer_type}*.bed")

#lFilelists=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/BarcodeGroup/LAML-KR.bed"]
    for idx, sPathFile in enumerate(lFilelists):
        data = pr.read_bed(sPathFile).df
        donors = np.unique(data["Score"]).tolist()
        num_donors = len(donors)
        i = num_donors - decrement_by
        while i > 0:
            os.makedirs(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}", exist_ok=True)
            for j in range(0, 100):
                subsampled_donors = random.sample(donors, i)
                data_subsampled = data.loc[data["Score"].isin(subsampled_donors), :]
                data_subsampled = pr.PyRanges(data_subsampled)
                subsample_fp = f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}/subsample_S{j}.bed"
                data_subsampled.to_bed(subsample_fp)
                outPath=f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}"
                os.makedirs(f"{outPath}/Intersected_paz_Cancergroup", exist_ok=True)
                os.makedirs(f"{outPath}/IntersectedCount_paz_Cancergroup", exist_ok=True)
                # outPath = f"{outPath}/Intersected_paz_Cancergroup"
                interval_fp = f"{current_dir}/../Sorted_Interval_paz.bed"
                IntersectBed(sPathFile, outPath, interval_fp, j)
            i = i - decrement_by







