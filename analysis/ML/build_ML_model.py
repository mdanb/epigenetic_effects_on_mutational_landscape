#### Imports ####
from ML_utils import *
from config import *
from natsort import natsorted
import subprocess

def run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type_or_donor_id, scATAC_dir,
                                         scATAC_cell_number_filter, annotation_dir,
                                         tissues_to_consider, ML_model, seed,
                                         n_optuna_trials_prebackward_selection,
                                         n_optuna_trials_backward_selection,
                                         tss_filter=None):
    #### Filter data ####
    scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir)
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    if (not pd.isna(mutations_df).any().any()):
        # for compatibility
        mutations_df = add_na_ranges(mutations_df)
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)
    cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type_or_donor_id)

    os.makedirs(f"models/{ML_model}/{cancer_type_or_donor_id}/{scATAC_dir}", exist_ok=True)

    if (tissues_to_consider == "all"):
        backwards_elim_dir=f"models/{ML_model}/" \
                           f"{cancer_type_or_donor_id}/{scATAC_dir}/backwards_elimination_results"
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
                       n_optuna_trials_backward_selection)

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
                                  n_optuna_trials_backward_selection, iters_dont_skip):
    start, end = map(int, seed_range.split('-'))
    seed_range = range(start, end + 1)
    for seed in seed_range:
        print(f"Running model for seed {seed}")
        scATAC_sources = construct_scATAC_sources(datasets)
        mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                      histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                                      CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor)
        if (not per_donor):
            for cancer_type in cancer_types:
                if (tss_fragment_filter):
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
                    run_unclustered_data_analysis_helper(datasets, mutations_df, cancer_type, scATAC_dir,
                                                         scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                         ML_model, seed, n_optuna_trials_prebackward_selection,
                                                         n_optuna_trials_backward_selection)
        else:
            scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter,
                                              tss_fragment_filter, annotation_dir, seed)
            for idx, donor in enumerate(mutations_df.columns):
                if (idx in range(*donor_range)):
                    run_unclustered_data_analysis_helper(datasets, mutations_df, donor, scATAC_dir,
                                                         scATAC_cell_number_filter, annotation_dir, tissues_to_consider,
                                                         ML_model, seed, n_optuna_trials_prebackward_selection,
                                                         n_optuna_trials_backward_selection)

        cancer_types = ",".join(cancer_types)
        datasets = ",".join(datasets)
        iters_dont_skip = ",".join(iters_dont_skip)

        subprocess.call(["Rscript", "plot_top_features.R",
                         f"--cancer_types={cancer_types}",
                         f"--ML_model={ML_model}",
                         f"--datasets={datasets}",
                         f"--seed={seed}",
                         f"--cell_number_filter={scATAC_cell_number_filter}",
                         f"--annotation={annotation_dir}",
                         f"--iters_dont_skip={iters_dont_skip}"])

    return


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
                                  meso, RNA_subtyped, per_donor, donor_range, ML_model, seed_range,
                                  n_optuna_trials_prebackward_selection, n_optuna_trials_backward_selection,
                                  iters_dont_skip)





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
