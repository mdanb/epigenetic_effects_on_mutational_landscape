import pyreadr
import pandas as pd
from natsort import natsorted

# Load Data helpers
def load_scATAC(scATAC_path):
    scATAC_df = pyreadr.read_r(scATAC_path)
    scATAC_df = scATAC_df[None]
    scATAC_df = scATAC_df.T
    return scATAC_df

def load_per_donor_mutations(cancer_type):
    df = pd.read_csv(f"../../data/processed_data/per_patient_mutations/{cancer_type}_per_donor.csv",
                     index_col=0)
    return(df.loc[natsorted(df.index)])

def load_agg_mutations():
    df = pd.read_csv("../../data/processed_data/mut_count_data.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_woo_pcawg_mutations():
    df = pd.read_csv("../../data/processed_data/pcawg_agg_woo.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])
# def load_meso_mutations(meso_waddell_and_biphasic, meso_waddell_only, meso_waddell_and_broad_only,
#                         meso_waddell_biph_786_846):
def load_meso_mutations(waddell_sarc_biph,
                        waddell_sarc,
                        waddell_sarc_tsankov_sarc,
                        waddell_sarc_biph_tsankov_sarc_biph):
    if (waddell_sarc_biph):
        # (waddell sarcomatoid + waddell biphasic = 5 samples) + waddell epithelioid
        df = pd.read_csv("../../data/processed_data/meso_mut_count_data.csv",
                           index_col=0)
    elif (waddell_sarc):
        # (just waddell sarcomatoid = 2 samples) + waddell epithelioid
        df = pd.read_csv("../../data/processed_data/mesothelioma_epithelioid_sarcomatoid2_WGS_Waddell.csv",
                           index_col=0)
    elif (waddell_sarc_tsankov_sarc):
        # (waddell sarcomatoid + tsankov sarcomatoid = 3 samples)
        df = pd.read_csv("../../data/processed_data/mesothelioma_p786_waddell_sarco.csv",
                           index_col=0)
    elif (waddell_sarc_biph_tsankov_sarc_biph):
         # (waddell sarcomatoid + waddell biphasic + tsankov sarcomatoid + tsankov biphasic
         # = 7 samples) + waddell epithelioid
         df = pd.read_csv("../../data/processed_data/mesothelioma_epithelioid_sarcomatoid_biphasic_WGS_Waddell_786_846.csv",
                           index_col=0)
    return(df.loc[natsorted(df.index)])

def load_sclc_mutations():
    df = pd.read_csv("../../data/processed_data/sclc_count_overlaps.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_subtyped_lung_mutations():
    df = pd.read_csv("../../data/processed_data/lung_subtyped_mutation.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_histologically_subtyped_mutations():
    df = pd.read_csv("../../data/processed_data/histologically_subtyped_mutations.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_de_novo_seurat_clustered_cancers(cancer_types):
    cancer_type = cancer_types[0].split("x")[0]
    seurat_cluster_settings = cancer_types[0].split("x")[1]
    df = pd.read_csv(f"../../data/processed_data/de_novo_seurat_clustered_mutations/{cancer_type}/"
                     f"{seurat_cluster_settings}.csv", index_col=0)
    # df = pd.read_csv("../../data/processed_data/de_novo_clustered_panc_adeno_mutations.csv",
    #                    index_col=0)
    return(df.loc[natsorted(df.index)])

# Filter Data Helpers
def filter_agg_data(scATAC_df, mutations_df):
    # if (meso):
    #     idx_select = scATAC_df.index.isin(mutations_df.index)
    #     scATAC_df = scATAC_df[idx_select]
    # else:
    idx_select = ~pd.isnull(mutations_df).any(axis=1)
    scATAC_df = scATAC_df.loc[idx_select]
    mutations_df = mutations_df[idx_select]
    return scATAC_df, mutations_df

def filter_mutations_by_cancer(mutations_df, cancer_type):
    cancer_type_specific_mut_idx = [i for i, s in enumerate(mutations_df.columns.values) if cancer_type == s][0]
    cancer_specific_mutations = mutations_df.iloc[:, cancer_type_specific_mut_idx]
    return cancer_specific_mutations

def add_na_ranges(mutations_df):
    full_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
    mutations_df = pd.merge(full_ranges, mutations_df, left_on = "x", right_index=True,
                            how="outer").set_index("x")
    return(mutations_df.loc[natsorted(mutations_df.index)])