#### Imports ####
from ML_utils import *
from config import *
from natsort import natsorted
import subprocess

# bioRxiv_method = config.bioRxiv_method
# def run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec, cancer_types,
#                                   meso_waddell_and_biphasic, meso_waddell_only, meso_waddell_and_broad_only,
#                                   meso_waddell_biph_786_846, tss_fragment_filter, scATAC_source="Bingren"):


def run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type_or_donor_id, scATAC_dir,
                                         scATAC_cell_number_filter, annotation_dir,
                                         tissues_to_consider, ML_model, seed, tss_filter=None):
    #### Filter data ####
    scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir)
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    if (not pd.isna(mutations_df).any().any()):
        # for compatibility
        mutations_df = add_na_ranges(mutations_df)
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)
    cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type_or_donor_id)

    os.makedirs(f"models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}",
                exist_ok=True)

    if (tissues_to_consider == "all"):
        backwards_elim_dir=f"models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results"
        grid_search_filename = f"models/{ML_model}/" \
                               f"{cancer_type_or_donor_id}/{scATAC_dir}/grid_search_results.pkl"
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance.txt"

        # All Cells
        train_val_test(scATAC_df, cancer_specific_mutations,
                       grid_search_filename,
                       backwards_elim_dir,
                       test_set_perf_filepath,
                       ML_model,
                       seed)

    # Tissue Specific
    else:
        tissues_string = "_".join(tissues_to_consider)
        backwards_elim_dir=f"models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results_{tissues_string}"
        grid_search_filename = f"models/{ML_model}/" \
                               f"{cancer_type_or_donor_id}/{scATAC_dir}/grid_search_results_{tissues_string}.pkl"
        test_set_perf_filepath = f"models/{ML_model}/" \
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
                       test_set_perf_filepath,
                       ML_model,
                       seed)

def run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                                  tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC,
                                  combined_CPTAC_ICGC, meso, RNA_subtyped, per_donor, donor_range, ML_model, seed):
    # waddell_sarc_biph_waddell_epith = config.waddell_sarc_biph_waddell_epith
    # waddell_sarc_waddell_epith = config.waddell_sarc_waddell_epith
    # waddell_sarc_tsankov_sarc_waddell_epith = config.waddell_sarc_tsankov_sarc_waddell_epith
    # waddell_sarc_biph_tsankov_sarc_biph = config.waddell_sarc_biph_tsankov_sarc_biph
    scATAC_sources = construct_scATAC_sources(datasets)
    # scATAC_dir_orig = f"scATAC_source_{scATAC_sources}_cell_number_filter_{scATAC_cell_number_filter}"
    mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                                  CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor)

    if (not per_donor):
        for cancer_type in cancer_types:
            if (tss_fragment_filter):
                for tss_filter in tss_fragment_filter:
                    scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_fragment_filter,
                                                      annotation_dir, seed)
                    # scATAC_dir = scATAC_dir_orig + "_tss_fragment_filter_" + tss_filter
                    run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir,
                                                         scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                         ML_model, seed, tss_filter=tss_filter)

            else:
                scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
                                                  tss_fragment_filter, annotation_dir, seed)
                run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir,
                                                     scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                     ML_model, seed)
    else:
        scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
                                          tss_fragment_filter, annotation_dir, seed)
        for idx, donor in enumerate(mutations_df.columns):
            if (idx in range(*donor_range)):
                run_unclustered_data_analysis_helper(datasets, mutations_df, donor, scATAC_dir,
                                                     scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                     ML_model, seed)
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
#         os.makedirs(f"analysis/ML/models/{cancer_type}", exist_ok=True)
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
#         cluster_dir = f"analysis/ML/models/{cancer_type}/{cluster_method_dir}/{threshold_dir}/{cluster_dir}/"
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

if (test_backward_selection_iter):
    save_n_features_model_test_performance(test_backward_selection_iter, datasets, ML_model, scATAC_cell_number_filter,
                                           tss_fragment_filter,
                                           annotation_dir, waddell_sarc_biph, waddell_sarc, waddell_sarc_tsankov_sarc,
                                           waddell_sarc_biph_tsankov_sarc_biph, SCLC, lung_subtyped, woo_pcawg,
                                           histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                                           CPTAC, combined_CPTAC_ICGC, per_donor, seed)
else:
    run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                                  tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC, combined_CPTAC_ICGC,
                                  meso, RNA_subtyped, per_donor, donor_range, ML_model, seed)
    cancer_types = ",".join(cancer_types)
    datasets = ",".join(datasets)
    subprocess.call(["Rscript", "plot_top_features.R",
                     f"--cancer_types={cancer_types}",
                     f"--ML_model={ML_model}",
                     f"--datasets={datasets}",
                     f"--seed={seed}",
                     f"--cell_number_filter={scATAC_cell_number_filter}",
                     f"--annotation={annotation_dir}"])





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
