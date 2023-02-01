import pandas as pd
import os
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cancer_types", nargs="+", type=str,
                    help='which cancer types to create mutation files for', required=True)
config = parser.parse_args()
cancer_types = config.cancer_types

def load_filenames_per_cancer(cancer_type):
    f = open(f"raw_dir/mutation_data/per_cancer_icgc/{cancer_type}/{cancer_type}.json")
    data = json.load(f)
    return data["files"]

def all_patients_are_processed(filepath, cancer_type):
    patients = pd.unique(pd.read_csv(filepath,
                         sep="\t", usecols=["icgc_donor_id"],
                         dtype={"icgc_donor_id":str})["icgc_donor_id"])
    for patient in patients:
        if (not os.path.exists(f"processed_data/per_patient_mutations/{cancer_type}/{patient}.csv")):
            return False
    return True

def load_and_clean_mutations(file):
    per_patient_mutations = pd.read_csv(file, sep="\t", usecols=["icgc_mutation_id", "icgc_donor_id",
                                                                 "chromosome", "chromosome_start",
                                                                 "chromosome_end"])
    unique_mutations = per_patient_mutations[["icgc_donor_id", "icgc_mutation_id"]].drop_duplicates()

    chromosome_start_end = per_patient_mutations[["icgc_mutation_id",
                                                   "icgc_donor_id",
                                                   "chromosome_start",
                                                   "chromosome_end",
                                                   "chromosome"]]
    chromosome_start_end = chromosome_start_end.drop_duplicates()

    unique_mutations = unique_mutations.merge(chromosome_start_end,
                           how="left",
                           on=["icgc_donor_id", "icgc_mutation_id"])
    return unique_mutations

def create_mutations_per_patient_files(mutations, cancer_type):
    os.makedirs(f"processed_data/per_patient_mutations/{cancer_type}", exist_ok=True)

    for index, group in mutations.groupby("icgc_donor_id"):
        group.to_csv(f"processed_data/per_patient_mutations/{cancer_type}/{index}.csv",
                    index=False)


# per_patient_mutations_lung_adeno = load_and_clean_mutations("raw_dir/mutation_data/per_cancer_icgc/LUAD-US/"
#                                                             "simple_somatic_mutation.open.LUAD-US.tsv")
# create_mutations_per_patient_files(per_patient_mutations_lung_adeno, "Lung-AdenoCA")

for cancer_type in cancer_types:
    files = load_filenames_per_cancer(cancer_type)
    for file in files:
        if not all_patients_are_processed(file, cancer_type):
            per_patient_mutations = load_and_clean_mutations(file)
            create_mutations_per_patient_files(per_patient_mutations, cancer_type)

# per_patient_mutations_skin_melanoma = load_and_clean_mutations("raw_dir/mutation_data/per_cancer_icgc/Skin-Melanoma/"
#                                                                 "simple_somatic_mutation.open.MELA-AU.tsv",
#                                                                "raw_dir/mutation_data/per_cancer_icgc/Skin-Melanoma/"
#                                                                 "simple_somatic_mutation.open.SKCM-US.tsv")