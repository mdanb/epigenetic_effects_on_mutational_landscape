#### Imports ####
from ML_utils import *
from config import create_parser
from natsort import natsorted
import glob

def run_unclustered_data_analysis_helper(scATAC_df, cancer_specific_mutations,
                                         cancer_type_or_donor_id, scATAC_dir,
                                         tissues_to_consider, ML_model, seed,
                                         n_optuna_trials_prebackward_selection,
                                         n_optuna_trials_backward_selection,
                                         backwards_elim_dir,
                                         feature_importance_method):

    os.makedirs(f"models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}", exist_ok=True)

    if tissues_to_consider == "all":
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance.txt"

        # All Cells
        train_val_test(scATAC_df, cancer_specific_mutations,
                       backwards_elim_dir,
                       test_set_perf_filepath,
                       ML_model,
                       seed,
                       scATAC_dir,
                       cancer_type_or_donor_id,
                       n_optuna_trials_prebackward_selection,
                       n_optuna_trials_backward_selection,
                       feature_importance_method)

    # Tissue Specific
    else:
        tissues_string = "_".join(tissues_to_consider)
        backwards_elim_dir=f"models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results_{tissues_string}"
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type_or_donor_id}/{scATAC_dir}/test_set_performance_{tissues_string}.txt"

        # tissue = cancer_type.split("-")[0]
        def check_tissue(tissue_cell_type):
            return(any(tissue in tissue_cell_type for tissue in tissues_to_consider))

        idx = pd.Series(scATAC_df.columns.values).apply(check_tissue)
        tissue_specific_cell_types = scATAC_df.columns.values[idx]
        per_tissue_df = scATAC_df.loc[:, tissue_specific_cell_types]
        train_val_test(per_tissue_df, cancer_specific_mutations,
                       backwards_elim_dir,
                       test_set_perf_filepath,
                       ML_model,
                       seed,
                       scATAC_dir,
                       cancer_type_or_donor_id,
                       n_optuna_trials_prebackward_selection,
                       n_optuna_trials_backward_selection)

def run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                                  tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC,
                                  combined_CPTAC_ICGC, meso, RNA_subtyped, per_donor, donor_range, ML_model,
                                  seed_range, n_optuna_trials_prebackward_selection,
                                  n_optuna_trials_backward_selection, top_features_to_plot, save_test_set_perf,
                                  make_plots, feature_importance_method):
    ### args used at the end for plot_top_features.R ###
    cancer_types_arg = ",".join(cancer_types)
    datasets_arg = ",".join(datasets)
    # iters_dont_skip_arg = ",".join(iters_dont_skip)
    ####################################################

    # mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
    #                               histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
    #                               CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor)
    scATAC_sources = construct_scATAC_sources(datasets)

    for seed in seed_range:
        print(f"Running model for seed {seed}")
        print(f"Using scATAC sources: {scATAC_sources}")
        if not per_donor:
            for cancer_type in cancer_types:
                if tss_fragment_filter:
                    for tss_filter in tss_fragment_filter:
                        scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_fragment_filter,
                                                          annotation_dir, seed)
                        run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir,
                                                             scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                             ML_model, seed, n_optuna_trials_prebackward_selection,
                                                             n_optuna_trials_backward_selection,
                                                             tss_filter=tss_filter)

                else:
                    scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
                                                      tss_fragment_filter, annotation_dir, seed)
                    print(f"scATAC_dir is {scATAC_dir}")
                    backwards_elim_dir=f"models/{ML_model}/" \
                    f"{cancer_type}/{scATAC_dir}/backwards_elimination_results"

                    scATAC_df, cancer_specific_mutations = load_data(meso, SCLC, lung_subtyped, woo_pcawg,
                                                                  histologically_subtyped_mutations,
                                                                  de_novo_seurat_clustering,
                                                                  CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
                                                                  datasets, scATAC_cell_number_filter,
                                                                  annotation_dir, cancer_type)

                    run_unclustered_data_analysis_helper(scATAC_df, cancer_specific_mutations,
                                                         cancer_type, scATAC_dir,
                                                         tissues_to_consider,
                                                         ML_model, seed, n_optuna_trials_prebackward_selection,
                                                         n_optuna_trials_backward_selection, backwards_elim_dir,
                                                         feature_importance_method)
                bp_path = f"../../figures/models/{ML_model}/{cancer_type}/{scATAC_dir}/" \
                          f"backwards_elimination_results/bar_plot.png"
                print(f"Bar plot path: {bp_path}")
            if make_plots:
                print("Making plots")
            # if not os.path.exists(bp_path):
                call_plot_top_features(seed, cancer_types_arg, ML_model, datasets_arg, scATAC_cell_number_filter,
                                       annotation_dir, top_features_to_plot, meso, SCLC, lung_subtyped, woo_pcawg,
                                       histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC,
                                       combined_CPTAC_ICGC, RNA_subtyped, per_donor)

            if save_test_set_perf:
                total_num_features = len(natsorted(glob.glob(f"{backwards_elim_dir}/*pkl"))) + 1
                for iter in range(total_num_features - 1):
                    if total_num_features - iter in str(top_features_to_plot):
                        save_iter_i_model_test_performance(iter + 1,
                                                           datasets, ML_model, scATAC_cell_number_filter,
                                                           tss_fragment_filter,
                                                           annotation_dir, meso, SCLC, lung_subtyped, woo_pcawg,
                                                           histologically_subtyped_mutations, de_novo_seurat_clustering,
                                                           cancer_types, CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
                                                           int(seed))

        else:
            scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
                                              tss_fragment_filter, annotation_dir, seed)
            for idx, donor in enumerate(mutations_df.columns):
                if idx in range(*donor_range):
                    run_unclustered_data_analysis_helper(datasets, mutations_df, donor, scATAC_dir,
                                                         scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                         ML_model, seed, n_optuna_trials_prebackward_selection,
                                                         n_optuna_trials_backward_selection)
    return


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
test_backward_selection_iters = config.test_backward_selection_iters
seed_range = config.seed_range
start, end = map(int, seed_range.split('-'))
seed_range = range(start, end + 1)
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
top_features_to_plot = config.top_features_to_plot
save_test_set_perf = config.save_test_set_perf
make_plots = config.make_plots
feature_importance_method = config.feature_importance_method

run_unclustered_data_analysis(datasets, cancer_types, scATAC_cell_number_filter, annotation_dir,
                               tissues_to_consider, tss_fragment_filter, SCLC, lung_subtyped, woo_pcawg,
                               histologically_subtyped_mutations, de_novo_seurat_clustering, CPTAC, combined_CPTAC_ICGC,
                               meso, RNA_subtyped, per_donor, donor_range, ML_model, seed_range,
                               n_optuna_trials_prebackward_selection, n_optuna_trials_backward_selection,
                               top_features_to_plot, save_test_set_perf, make_plots, feature_importance_method)





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
