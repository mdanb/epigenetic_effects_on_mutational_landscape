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

def IntersectBed(cancer_type, outPath, interval_fp, subsample_idx):
    sOutFile = f"{cancer_type}_S{subsample_idx}.bed"
    sInFile = f"{outPath}/subsample_S{subsample_idx}.bed"
    os.system(f"bedtools intersect -wa -a {interval_fp} -b " + sInFile + f" -loj > {outPath}/Intersected_paz_Cancergroup/Intersected_"+sOutFile)
   # sleep(1)
    os.system(f"bedtools intersect -wa -c -a {interval_fp} -b " + sInFile + f" -loj > {outPath}/IntersectedCount_paz_Cancergroup/IntersectedCount_"+sOutFile)
    #sleep(1)

if __name__=="__main__":
    lFilelists = []
    for cancer_type in cancer_types:
        os.makedirs(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/", exist_ok=True)
        lFilelists = lFilelists + glob.glob(f"{current_dir}/../mutation_data/bed_files/{cancer_type}*.bed")
    print(lFilelists)
#lFilelists=["/ahg/regevdata/projects/ICA_Lung/Wooseung/CellOrigin/Data/BarcodeGroup/LAML-KR.bed"]
    for idx, sPathFile in enumerate(lFilelists):
        data = pr.read_bed(sPathFile).df
        donors = np.unique(data["Score"]).tolist()
        num_donors = len(donors)
        i = num_donors - decrement_by
        while i > 0:
            os.makedirs(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}", exist_ok=True)
            for j in range(0, 100):
                subsample_fp = f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}/subsample_S{j}.bed"
                if not os.path.exists(f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}/subsample_S{j+1}.bed"):
                    subsampled_donors = random.sample(donors, i)
                    data_subsampled = data.loc[data["Score"].isin(subsampled_donors), :]
                    data_subsampled = pr.PyRanges(data_subsampled)
                    data_subsampled.to_bed(subsample_fp)
                    outPath=f"{current_dir}/../mutation_data/bed_files/{cancer_type}_subsampled/{i}"
                    os.makedirs(f"{outPath}/Intersected_paz_Cancergroup", exist_ok=True)
                    os.makedirs(f"{outPath}/IntersectedCount_paz_Cancergroup", exist_ok=True)
                    # outPath = f"{outPath}/Intersected_paz_Cancergroup"
                    interval_fp = f"{current_dir}/../Sorted_Interval_paz.bed"
                    IntersectBed(cancer_type, outPath, interval_fp, j)
            i = i - decrement_by







