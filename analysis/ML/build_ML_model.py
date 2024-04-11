#### Imports ####

from ML_utils import *
from config import create_parser
import time

def run_unclustered_data_analysis_helper(scATAC_df,
                                         cancer_specific_mutations,
                                         cancer_type_or_donor_id,
                                         scATAC_sources,
                                         ML_model,
                                         seed,
                                         n_optuna_trials_prebackward_selection,
                                         n_optuna_trials_backward_selection,
                                         feature_importance_method,
                                         sqlite,
                                         debug_bfs,
                                         fold_for_test_set,
                                         hundred_kb,
                                         expanded_hundred_kb,
                                         scATAC_cell_number_filter,
                                         annotation_dir,
                                         tss_fragment_filter,
                                         # save_test_set_perf,
                                         test_set_perf_num_features,
                                         tissues_to_consider,
                                         grid_analysis,
                                         grid_cell_type,
                                         cell_types_keep):

    tissues_string = "_".join(tissues_to_consider)
    # if cell_types_keep:
    #     cell_types_keep = "_".join(cell_types_keep)
    scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_fragment_filter, annotation_dir,
                                      hundred_kb, expanded_hundred_kb, tissues_string, seed=seed,
                                      fold_for_test_set=fold_for_test_set, grid_analysis=grid_analysis,
                                      cell_types_keep=cell_types_keep)

    print(grid_cell_type)
    if grid_analysis:
        if grid_cell_type[0] == "cerebrum Astrocytes/Oligodendrocytes SH":
           grid_cell_type[0] = "cerebrum Astrocytes-Oligodendrocytes SH"
        backwards_elim_dir=f"models/{ML_model}/{cancer_type_or_donor_id}_{grid_cell_type[0].replace(' ', '-')}/" \
                           f"{scATAC_dir}/backwards_elimination_results"
    else:
        backwards_elim_dir=f"models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results"

    backwards_elim_dir = f"{os.path.dirname(os.path.realpath(__file__))}/{backwards_elim_dir}"
    os.makedirs(f"{backwards_elim_dir}", exist_ok=True)

    # if tissues_to_consider == "all":
    if grid_analysis:
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}_{grid_cell_type[0].replace(' ', '-')}/{scATAC_dir}/" \
                                 f"test_set_performance.txt"
        cancer_type_or_donor_id = f"{cancer_type_or_donor_id}_{grid_cell_type[0].replace(' ', '-')}"
    else:
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance.txt"

    start = time.time()
    train_val_test(scATAC_df,
                   cancer_specific_mutations,
                   backwards_elim_dir,
                   test_set_perf_filepath,
                   ML_model,
                   seed,
                   scATAC_dir,
                   cancer_type_or_donor_id,
                   n_optuna_trials_prebackward_selection,
                   n_optuna_trials_backward_selection,
                   feature_importance_method,
                   sqlite,
                   debug_bfs,
                   fold_for_test_set,
                   hundred_kb)
    end = time.time()
    print(f"{(end - start) / 60} minutes")
    time_fp = f"{os.path.dirname(os.path.realpath(__file__))}/models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}/time.txt"
    with open(time_fp, "a") as f:
        f.write(f"{(end - start) / 60}\n")
    print(f"Done modeling {cancer_type_or_donor_id}!")

    # if save_test_set_perf:
    print("Saving test set performances for backward elimination...")
    # total_num_features = len(natsorted(glob.glob(f"{backwards_elim_dir}/*pkl"))) + 1
    # for curr_num_feats in range(1, total_num_features):
        # if curr_num_feats in test_set_perf_num_features:
    #range(1, total_num_features)
    if test_set_perf_num_features[0] == "all":
        if scATAC_df.shape[1] > 20:
            test_set_perf_num_features = [i for i in range(1, 20)]
        else:
            test_set_perf_num_features = [i for i in range(1, scATAC_df.shape[1] + 1)]

    print(test_set_perf_num_features)
    for curr_num_feats in test_set_perf_num_features:
       save_model_with_n_features_test_performance(scATAC_df, cancer_specific_mutations, scATAC_dir, curr_num_feats,
                                                   ML_model, cancer_type_or_donor_id, feature_importance_method,
                                                   fold_for_test_set)
    print("Done saving test set performances!")

    # # Tissue Specific
    # else:
    #     tissues_string = "_".join(tissues_to_consider)
    #     backwards_elim_dir=f"models/{ML_model}/" \
    #                        f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results_{tissues_string}"
    #     test_set_perf_filepath = f"models/{ML_model}/" \
    #                              f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance_{tissues_string}.txt"
    #
    #     # tissue = cancer_type.split("-")[0]
    #     def check_tissue(tissue_cell_type):
    #         return any(tissue in tissue_cell_type for tissue in tissues_to_consider)
    #
    #     idx = pd.Series(scATAC_df.columns.values).apply(check_tissue)
    #     tissue_specific_cell_types = scATAC_df.columns.values[idx]
    #     per_tissue_df = scATAC_df.loc[:, tissue_specific_cell_types]
    #     train_val_test(per_tissue_df, cancer_specific_mutations,
    #                    backwards_elim_dir,
    #                    test_set_perf_filepath,
    #                    ML_model,
    #                    seed,
    #                    scATAC_dir,
    #                    cancer_type_or_donor_id,
    #                    n_optuna_trials_prebackward_selection,
    #                    n_optuna_trials_backward_selection)

def run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                                  tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC,
                                  combined_CPTAC_ICGC, meso, RNA_subtyped, hierarchically_subtyped_mutations, mm,
                                  msi_high, hundred_kb, expanded_hundred_kb,
                                  per_donor, donor_range, aggregated_per_donor, ML_model, seed_range,
                                  n_optuna_trials_prebackward_selection,
                                  n_optuna_trials_backward_selection, top_features_to_plot,
                                  make_plots, feature_importance_method, sqlite, test_set_perf_num_features,
                                  debug_bfs, fold_for_test_set, tissues_to_consider, grid_analysis, grid_cell_types,
                                  cell_types_keep, subsampled_mutations, custom_mutations, which_interval_ranges,
                                  dataset_abbrev):
    ### args used at the end for plot_top_features.R ###
    scATAC_sources = construct_scATAC_sources(datasets)
    scATAC_df = construct_scATAC_df(tss_fragment_filter, datasets, scATAC_cell_number_filter, annotation_dir,
                                    hundred_kb, expanded_hundred_kb, tissues_to_consider,
                                    grid_analysis, grid_cell_types, cell_types_keep, which_interval_ranges,
                                    dataset_abbrev)
    if not grid_analysis:
        grouped_cell_types = [scATAC_df.columns]
    else:
        # 1 cell type per group i.e cell type grid analysis
        grouped_cell_types = grid_cell_types

    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    # Note that the loading process arranges the bins so that it's in the correct order of the genome
    # ensuring that splits truly split based on contiguous genomic regions
    for seed in seed_range:
        print(f"Running model for seed {seed}")
        print(f"Using scATAC sources: {scATAC_sources}")
        for cancer_type in cancer_types:
            for cell_type_group in grouped_cell_types:
                if isinstance(cell_type_group, str):
                    cell_type_group = [cell_type_group]
                scdf = scATAC_df.loc[:, cell_type_group]
                print(f"Working on {cancer_type}...")
                mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                              histologically_subtyped_mutations, de_novo_seurat_clustering,
                                              CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, cancer_type,
                                              hundred_kb, expanded_hundred_kb, aggregated_per_donor,
                                              hierarchically_subtyped_mutations, mm, msi_high, subsampled_mutations,
                                              custom_mutations)
                if subsampled_mutations:
                    cancer_type_name_for_subsample = "_".join([cancer_type, "seed", str(seed), "fold",
                                                               str(fold_for_test_set + 1)])
                else:
                    cancer_type_name_for_subsample = None

                scdf, cancer_specific_mutations = prep_and_align_mutations_with_scatac(scdf,
                                                                                        mutations_df,
                                                                                        cancer_type,
                                                                                        hundred_kb,
                                                                                        expanded_hundred_kb,
                                                                                        per_donor,
                                                                                        subsampled_mutations,
                                                                                        cancer_type_name_for_subsample)
                if not per_donor:
                    run_unclustered_data_analysis_helper(scdf,
                                                         cancer_specific_mutations,
                                                         cancer_type,
                                                         scATAC_sources,
                                                         ML_model,
                                                         seed,
                                                         n_optuna_trials_prebackward_selection,
                                                         n_optuna_trials_backward_selection,
                                                         feature_importance_method,
                                                         sqlite,
                                                         debug_bfs,
                                                         fold_for_test_set,
                                                         hundred_kb,
                                                         expanded_hundred_kb,
                                                         scATAC_cell_number_filter,
                                                         annotation_dir,
                                                         tss_fragment_filter,
                                                         # save_test_set_perf,
                                                         test_set_perf_num_features,
                                                         tissues_to_consider,
                                                         grid_analysis,
                                                         cell_type_group,
                                                         cell_types_keep)
                else:
                    for idx, donor in enumerate(mutations_df.columns):
                        if idx in range(*donor_range):
                            donor_specific_mutations = mutations_df[donor]
                            run_unclustered_data_analysis_helper(scdf,
                                                                 donor_specific_mutations,
                                                                 f"{cancer_type}_{donor}",
                                                                 scATAC_sources,
                                                                 ML_model,
                                                                 seed,
                                                                 n_optuna_trials_prebackward_selection,
                                                                 n_optuna_trials_backward_selection,
                                                                 feature_importance_method,
                                                                 sqlite,
                                                                 debug_bfs,
                                                                 fold_for_test_set,
                                                                 hundred_kb,
                                                                 expanded_hundred_kb,
                                                                 scATAC_cell_number_filter,
                                                                 annotation_dir,
                                                                 tss_fragment_filter,
                                                                 # save_test_set_perf,
                                                                 test_set_perf_num_features,
                                                                 tissues_to_consider)

        if make_plots:
            call_plot_top_features(seed_range, cancer_types, ML_model, datasets, scATAC_cell_number_filter,
                                   annotation_dir, top_features_to_plot, feature_importance_method, tissues_to_consider,
                                   fold_for_test_set + 1)
    return


parser = create_parser()
config = parser.parse_args()
cancer_types = config.cancer_types
# run_clustered_mutations = config.clustered_mutations
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
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
meso = config.meso
RNA_subtyped = config.RNA_subtyped
per_donor = config.per_donor
donor_range = config.donor_range
if per_donor:
    assert not (donor_range is None)
seed_range = config.seed_range
start, end = map(int, seed_range.split('-'))
seed_range = range(start, end + 1)
fold_for_test_set = config.fold_for_test_set
assert 1 <= fold_for_test_set <= 10
fold_for_test_set = fold_for_test_set - 1
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
top_features_to_plot = config.top_features_to_plot
# save_test_set_perf = config.save_test_set_perf
make_plots = config.make_plots
feature_importance_method = config.feature_importance_method
sqlite = config.sqlite
test_set_perf_num_features = config.test_set_perf_num_features
debug_bfs = config.debug_bfs
hundred_kb = config.hundred_kb
expanded_hundred_kb = config.expanded_hundred_kb
aggregated_per_donor = config.aggregated_per_donor
hierarchically_subtyped_mutations = config.hierarchically_subtyped_mutations
grid_analysis = config.grid_analysis
mm = config.mm
msi_high = config.msi_high
cell_types_keep = config.cell_types_keep
grid_cell_types = None
if grid_analysis:
    grid_cell_types = config.grid_cell_types.split(",")
subsampled_mutations = config.subsampled_mutations
custom_mutations = config.custom_mutations
which_interval_ranges = config.which_interval_ranges
dataset_abbrev = config.dataset_abbrev

run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                               tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                               histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC, combined_CPTAC_ICGC,
                               meso, RNA_subtyped, hierarchically_subtyped_mutations, mm, msi_high,
                               hundred_kb, expanded_hundred_kb, per_donor, donor_range, aggregated_per_donor, ML_model,
                               seed_range, n_optuna_trials_prebackward_selection, n_optuna_trials_backward_selection,
                               top_features_to_plot, make_plots, feature_importance_method,
                               sqlite, test_set_perf_num_features, debug_bfs, fold_for_test_set, tissues_to_consider,
                               grid_analysis, grid_cell_types, cell_types_keep, subsampled_mutations, custom_mutations,
                               which_interval_ranges, dataset_abbrev)



# scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
#                                           tss_fragment_filter,
#                                           annotation_dir, hundred_kb,
#                                           seed, fold_for_test_set)


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


# scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
#                                   tss_fragment_filter,
#                                   annotation_dir, hundred_kb,
#                                   seed, fold_for_test_set)
# backwards_elim_dir=f"models/{ML_model}/" \
# f"{cancer_type}/{scATAC_dir}/backwards_elimination_results"

# # Note that the loading process arranges the bins so that it's in the correct order of the genome
# # ensuring that splits truly split based on contiguous genomic regions
# mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
#                               histologically_subtyped_mutations, de_novo_seurat_clustering,
#                               CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, cancer_type,
#                               hundred_kb)
#
# scATAC_df, cancer_specific_mutations = prep_and_align_mutations_with_scatac(scATAC_df, mutations_df,
#                                                                             cancer_type,
#                                                                             hundred_kb)
# scATAC_df, cancer_specific_mutations = load_data(meso, SCLC, lung_subtyped, woo_pcawg,
#                                                   histologically_subtyped_mutations,
#                                                   de_novo_seurat_clustering,
#                                                   CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
#                                                   datasets, scATAC_cell_number_filter,
#                                                   annotation_dir, cancer_type, hundred_kb)


# if test_backward_selection_iter:
#     for seed in seed_range:
#         try:
#             # Not all runs may be done
#             save_n_features_model_test_performance(test_backward_selection_iter, datasets, ML_model, scATAC_cell_number_filter,
#                                                    tss_fragment_filter,
#                                                    annotation_dir, meso, SCLC, lung_subtyped, woo_pcawg,
#                                                    histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
#                                                    CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, int(seed))
#         except:
#             pass
# else:


# run_unclustered_data_analysis_helper(scATAC_df, cancer_specific_mutations,
#                                      cancer_type, scATAC_dir,
#                                      tissues_to_consider,
#                                      ML_model, seed, n_optuna_trials_prebackward_selection,
#                                      n_optuna_trials_backward_selection, backwards_elim_dir,
#                                      feature_importance_method, sqlite, debug_bfs,
#                                      fold_for_test_set, hundred_kb)

# run_unclustered_data_analysis_helper(datasets, mutations_df, donor, scATAC_dir,
#                                      scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
#                                      ML_model, seed, n_optuna_trials_prebackward_selection,
#                                      n_optuna_trials_backward_selection)


# cancer_types_arg = ",".join(cancer_types)
# datasets_arg = ",".join(datasets)


# if make_plots:
# # if not os.path.exists(bp_path):
#     call_plot_top_features(seed, cancer_types_arg, ML_model, datasets_arg, scATAC_cell_number_filter,
#                            annotation_dir, top_features_to_plot, feature_importance_method,
#                            fold_for_test_set)
