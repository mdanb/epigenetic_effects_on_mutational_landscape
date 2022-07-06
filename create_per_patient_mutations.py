import pandas as pd
import os

def load_and_clean_mutations(filepath, *args):
    per_patient_mutations = pd.read_csv(filepath, sep="\t")
    for file in args:
        df = pd.read_csv(file, sep="\t")
        per_patient_mutations = pd.concat((per_patient_mutations, df), axis=0)
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


per_patient_mutations_lung_adeno = load_and_clean_mutations("raw_dir/mutation_data/per_cancer_icgc/LUAD-US/"
                                                            "simple_somatic_mutation.open.LUAD-US.tsv")
per_patient_mutations_skin_melanoma = load_and_clean_mutations("raw_dir/mutation_data/per_cancer_icgc/Skin-Melanoma/"
                                                                "simple_somatic_mutation.open.MELA-AU.tsv",
                                                               "raw_dir/mutation_data/per_cancer_icgc/Skin-Melanoma/"
                                                                "simple_somatic_mutation.open.SKCM-US.tsv")
create_mutations_per_patient_files(per_patient_mutations_lung_adeno, "Lung-AdenoCA")
create_mutations_per_patient_files(per_patient_mutations_skin_melanoma, "Skin-Melanoma")