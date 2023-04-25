#### Imports ####
import argparse
from ML_utils import *
from natsort import natsorted

def range_type(range_str):
    start, end = map(int, range_str.split('-'))
    if start > end:
        raise argparse.ArgumentTypeError('Invalid range: start value must be '
                                         'less than or equal to end value')
    return (start, end)

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', default=None)
parser.add_argument('--clustered_mutations', action="store_true",
                    help='run model on hierarchically clustered mutations', default=False)
parser.add_argument('--datasets', nargs="+", type=str,
                    help='which sc-ATACseq datasets to analyze', required=True)
parser.add_argument('--scATAC_cell_number_filter', type=int,
                    help='minimum number of cells per cell type in scATAC', default=100)
parser.add_argument('--annotation_dir', type=str,
                    help='name of annotation directory', default="default_annotation")
group = parser.add_mutually_exclusive_group()
group.add_argument('--waddell_sarc_biph', action="store_true",
                    default=False)
group.add_argument('--waddell_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_tsankov_sarc', action="store_true", default=False)
group.add_argument('--waddell_sarc_biph_tsankov_sarc_biph', action="store_true", default=False)
group.add_argument("--SCLC", action="store_true", default=False)
group.add_argument("--lung_subtyped", action="store_true", default=False)
parser.add_argument("--woo_pcawg", action="store_true", default=False)
parser.add_argument("--histologically_subtyped_mutations", action="store_true", default=False)
parser.add_argument("--de_novo_seurat_clustering", action="store_true", default=False)
parser.add_argument("--per_donor", action="store_true", default=False)
parser.add_argument("--CPTAC", action="store_true", default=False)
parser.add_argument("--combined_CPTAC_ICGC", action="store_true", default=False)
parser.add_argument('--donor_range', type=range_type, help='Specify a range in the format start-end',
                    default=None)
parser.add_argument('--tss_fragment_filter', nargs="+", type=str,
                    help='tss fragment filters to consider', default="")
parser.add_argument('--tissues_to_consider', nargs="+", type=str, default="all")
parser.add_argument("--ML_model", type=str, default="RF")

config = parser.parse_args()
cancer_types = config.cancer_types
run_clustered_mutations = config.clustered_mutations
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
waddell_sarc_biph = config.waddell_sarc_biph
waddell_sarc = config.waddell_sarc
waddell_sarc_tsankov_sarc = config.waddell_sarc_tsankov_sarc
waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
annotation_dir = config.annotation_dir
SCLC = config.SCLC
CPTAC = config.CPTAC
combined_CPTAC_ICGC = config.combined_CPTAC_ICGC
tissues_to_consider = config.tissues_to_consider
# tss_filtered = config.tss_filtered
tss_fragment_filter = config.tss_fragment_filter
ML_model = config.ML_model
lung_subtyped = config.lung_subtyped
woo_pcawg = config.woo_pcawg
histologically_subtyped_mutations = config.histologically_subtyped_mutations
de_novo_seurat_clustering = config.de_novo_seurat_clustering
per_donor = config.per_donor
donor_range = config.donor_range

# bioRxiv_method = config.bioRxiv_method
# def run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec, cancer_types,
#                                   meso_waddell_and_biphasic, meso_waddell_only, meso_waddell_and_broad_only,
#                                   meso_waddell_biph_786_846, tss_fragment_filter, scATAC_source="Bingren"):


def run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type_or_donor_id, scATAC_dir,
                                         waddell_sarc_biph, waddell_sarc, waddell_sarc_tsankov_sarc,
                                         waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter, annotation_dir,
                                         tissues_to_consider, ML_model, tss_filter=None):
    #### Filter data ####
    scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir)
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    if (not pd.isna(mutations_df).any().any()):
        # for compatibility
        mutations_df = add_na_ranges(mutations_df)
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)
    cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type_or_donor_id)

    scATAC_dir = append_meso_to_dirname_as_necessary(waddell_sarc_biph, waddell_sarc, waddell_sarc_tsankov_sarc,
                                                     waddell_sarc_biph_tsankov_sarc_biph, scATAC_dir)
    scATAC_dir = scATAC_dir + f"_annotation_{annotation_dir}"

    os.makedirs(f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}",
                exist_ok=True)

    if (tissues_to_consider == "all"):
        backwards_elim_dir=f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results"
        grid_search_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                               f"{cancer_type_or_donor_id}/{scATAC_dir}/grid_search_results.pkl"
        test_set_perf_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance.txt"

        # All Cells
        train_val_test(scATAC_df, cancer_specific_mutations,
                       grid_search_filename,
                       backwards_elim_dir,
                       test_set_perf_filename,
                       ML_model)

    # Tissue Specific
    else:
        tissues_string = "_".join(tissues_to_consider)
        backwards_elim_dir=f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results_{tissues_string}"
        grid_search_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                               f"{cancer_type_or_donor_id}/{scATAC_dir}/grid_search_results_{tissues_string}.pkl"
        test_set_perf_filename = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance_{tissues_string}.txt"

        # tissue = cancer_type.split("-")[0]
        def check_tissue(tissue_cell_type):
            return(any(tissue in tissue_cell_type for tissue in tissues_to_consider))

        idx = pd.Series(scATAC_df.columns.values).apply(check_tissue)
        tissue_specific_cell_types = scATAC_df.columns.values[idx]
        per_tissue_df = scATAC_df.loc[:, tissue_specific_cell_types]
        train_val_test(per_tissue_df,
                       cancer_specific_mutations,
                       grid_search_filename,
                       backwards_elim_dir,
                       test_set_perf_filename)

def run_unclustered_data_analysis(datasets, cancer_types, waddell_sarc_biph, waddell_sarc, waddell_sarc_tsankov_sarc,
                                  waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter, annotation_dir,
                                  tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC, combined_CPTAC_ICGC,
                                  per_donor, donor_range, ML_model):
    # waddell_sarc_biph_waddell_epith = config.waddell_sarc_biph_waddell_epith
    # waddell_sarc_waddell_epith = config.waddell_sarc_waddell_epith
    # waddell_sarc_tsankov_sarc_waddell_epith = config.waddell_sarc_tsankov_sarc_waddell_epith
    # waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
    scATAC_sources = ""

    for idx, dataset in enumerate(datasets):
        if (scATAC_sources == ""):
            scATAC_sources = dataset
        else:
            scATAC_sources = "_".join((scATAC_sources, dataset))

    scATAC_dir_orig = f"scATAC_source_{scATAC_sources}_cell_number_filter_{scATAC_cell_number_filter}"
    if (waddell_sarc_biph or waddell_sarc or waddell_sarc_tsankov_sarc or
        waddell_sarc_biph_tsankov_sarc_biph):
        mutations_df = load_meso_mutations(waddell_sarc_biph,
                                           waddell_sarc,
                                           waddell_sarc_tsankov_sarc,
                                           waddell_sarc_biph_tsankov_sarc_biph)
    elif (SCLC):
        mutations_df = load_sclc_mutations()
    elif (lung_subtyped):
        mutations_df = load_subtyped_lung_mutations()
    elif (woo_pcawg):
        mutations_df = load_woo_pcawg_mutations()
    elif (histologically_subtyped_mutations):
        mutations_df = load_histologically_subtyped_mutations()
    elif (de_novo_seurat_clustering):
        mutations_df = load_de_novo_seurat_clustered_cancers(cancer_types)
    elif (CPTAC):
        mutations_df = load_CPTAC()
    elif (combined_CPTAC_ICGC):
        mutations_df = load_combined_CPTAC_ICGC()
    elif (per_donor):
        mutations_df = load_per_donor_mutations(cancer_types[0])
    else:
        mutations_df = load_agg_mutations()

    if (not per_donor):
        for cancer_type in cancer_types:
            if (tss_fragment_filter):
                for tss_filter in tss_fragment_filter:
                    scATAC_dir = scATAC_dir_orig + "_tss_fragment_filter_" + tss_filter
                    run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir, waddell_sarc_biph,
                                                         waddell_sarc, waddell_sarc_tsankov_sarc,
                                                         waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter,
                                                         annotation_dir, tissues_to_consider, ML_model,
                                                         tss_filter=tss_filter)

            else:
                run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir_orig, waddell_sarc_biph,
                                                     waddell_sarc, waddell_sarc_tsankov_sarc,
                                                     waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter,
                                                     annotation_dir, tissues_to_consider, ML_model)
    else:
        for idx, donor in enumerate(mutations_df.columns):
            if (idx in range(*donor_range)):
                run_unclustered_data_analysis_helper(datasets, mutations_df, donor, scATAC_dir_orig,
                                                     waddell_sarc_biph,
                                                     waddell_sarc, waddell_sarc_tsankov_sarc,
                                                     waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter,
                                                     annotation_dir, tissues_to_consider, ML_model)
        # if (bioRxiv_method):
        #     tissues = set(scATAC_df.columns.str.split().to_series().apply(lambda x: x[0]))
        #     for tissue in tissues:
        #         grid_search_filename = f"models/bioRxiv_method/{cancer_type}/{scATAC_dir}/grid_search_results.pkl"
        #         test_set_perf_filename = f"models/bioRxiv_method/{cancer_type}/{scATAC_dir}/test_set_performance.txt"
        #
        #         train_val_test(scATAC_df, cancer_specific_mutations,
        #                        grid_search_filename,
        #                        backwards_elim_dir,
        #                        test_set_perf_filename)


# def run_clustered_data_analysis(scATAC_df, cancer_types):
#     for cancer_type in cancer_types:
#         os.makedirs(f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}", exist_ok=True)
#         cancer_hierarchical_dir = f"../../data/processed_data/hierarchically_clustered_mutations/{cancer_type}"
#         for cluster_method_dir in os.listdir(cancer_hierarchical_dir):
#             for threshold_dir in os.listdir(f"{cancer_hierarchical_dir}/{cluster_method_dir}"):
#                 run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
#                                        threshold_dir)
#
# def run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
#                            threshold_dir):
#     for cluster_file in glob.glob(f"{cancer_hierarchical_dir}/{cluster_method_dir}/{threshold_dir}/aggregated*"):
#         mutations_df = pd.read_csv(cluster_file, index_col=0)
#         cluster_dir = cluster_file.split("/")[-1].split(".")[0]
#         cluster_dir = f"/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ML/models/{cancer_type}/{cluster_method_dir}/{threshold_dir}/{cluster_dir}/"
#         os.makedirs(cluster_dir, exist_ok=True)
#         backwards_elim_dir=f"{cluster_dir}/backwards_elimination_results"
#         grid_search_filename = f"{cluster_dir}/grid_search_results.pkl"
#         test_set_perf_filename = f"{cluster_dir}/test_set_performance.txt"
#
#         scATAC_df = filter_clustered_data(scATAC_df, mutations_df)
#         mutations = mutations_df.values.reshape(-1)
#         train_val_test(scATAC_df, mutations,
#                        grid_search_filename,
#                        backwards_elim_dir,
#                        test_set_perf_filename)


# scATAC_df = pd.concat((scATAC_df, scATAC_df_bingren), axis=1)
# scATAC_sources = scATAC_sources + "bing_ren"
# scATAC_df = pd.concat((scATAC_df, scATAC_df_shendure), axis=1)
# scATAC_sources = scATAC_sources + "shendure"
# scATAC_df = pd.concat((scATAC_df, scATAC_df_tsankov), axis=1)
# scATAC_sources = scATAC_sources + "tsankov"

# if (len(datasets) == 5):
#     scATAC_sources = "combined_datasets"

# if (run_all_cells or run_tissue_spec_cells):

run_unclustered_data_analysis(datasets, cancer_types, waddell_sarc_biph, waddell_sarc, waddell_sarc_tsankov_sarc,
                              waddell_sarc_biph_tsankov_sarc_biph, scATAC_cell_number_filter, annotation_dir,
                              tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                              histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC, combined_CPTAC_ICGC,
                              per_donor, donor_range, ML_model)


    # run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec_cells,
    #                               cancer_types, waddell_sarc_biph, waddell_sarc,
    #                               waddell_sarc_tsankov_sarc, waddell_sarc_biph_tsankov_sarc_biph,
    #                               tss_fragment_filter, scATAC_sources)


# if (run_all_cells and shendure):
#     run_unclustered_data_analysis(scATAC_df_shendure, run_all_cells, run_tissue_spec_cells, cancer_types, meso,
#                                   "shendure")


# if (run_clustered_mutations and bing_ren):
#     run_clustered_data_analysis(scATAC_df, cancer_types)

# if (bioRxiv_method):
#     run_unclustered_data_analysis(combined_scATAC_df, run_all_cells, run_tissue_spec_cells,
#                                   cancer_types, meso, "combined_datasets")
