import argparse

def range_type(range_str):
    start, end = map(int, range_str.split('-'))
    if start > end:
        raise argparse.ArgumentTypeError('Invalid range: start value must be '
                                         'less than or equal to end value')
    return (start, end)

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cancer_types', nargs="+", type=str,
                        help='which cancer types to analyze', default=None)
    # parser.add_argument('--clustered_mutations', action="store_true",
    #                     help='run model on hierarchically clustered mutations', default=False)
    parser.add_argument('--datasets', nargs="+", type=str,
                        help='which sc-ATACseq datasets to analyze', required=True)
    parser.add_argument('--scATAC_cell_number_filter', type=int,
                        help='minimum number of cells per cell type in scATAC', default=100)
    parser.add_argument('--annotation_dir', type=str,
                        help='name of annotation directory', default="default_annotation")
    parser.add_argument('--tss_fragment_filter', nargs="+", type=str,
                        help='tss fragment filters to consider', default="")
    parser.add_argument('--tissues_to_consider', nargs="+", type=str, default="all")
    parser.add_argument("--ML_model", type=str, default="XGB")
    parser.add_argument('--test_backward_selection_iters', type=int, nargs="+", default=None)
    parser.add_argument('--seed_range', type=str)
    parser.add_argument('--n_optuna_trials_prebackward_selection', type=int, default=None)
    parser.add_argument('--n_optuna_trials_backward_selection', type=int, default=None)
    parser.add_argument('--top_features_to_plot', nargs="+", type=int)
    parser.add_argument("--save_test_set_perf", action="store_true", default=False)
    parser.add_argument("--make_plots", action="store_true", default=False)
    parser.add_argument("--feature_importance_method", type=str, default="permutation_importance")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--SCLC", action="store_true", default=False)
    group.add_argument("--lung_subtyped", action="store_true", default=False)
    group.add_argument("--woo_pcawg", action="store_true", default=False)
    group.add_argument("--histologically_subtyped_mutations", action="store_true", default=False)
    group.add_argument("--de_novo_seurat_clustering", action="store_true", default=False)
    group.add_argument("--per_donor", action="store_true", default=False)
    group.add_argument("--CPTAC", action="store_true", default=False)
    group.add_argument("--meso", action="store_true", default=False)
    group.add_argument("--combined_CPTAC_ICGC", action="store_true", default=False)
    group.add_argument('--donor_range', type=range_type, help='Specify a range in the format start-end',
                        default=None)
    group.add_argument("--RNA_subtyped", action="store_true", default=False)
    return parser
# parser.add_argument('--iters_dont_skip', nargs="+", type=str, default=["18"])
# parser.add_argument('--num_iter_skips', type=int, default=5)

# def extract_config(config):
#     cancer_types = config.cancer_types
#     datasets = sorted(config.datasets)
#     scATAC_cell_number_filter = config.scATAC_cell_number_filter
#     annotation_dir = config.annotation_dir
#     SCLC = config.SCLC
#     CPTAC = config.CPTAC
#     combined_CPTAC_ICGC = config.combined_CPTAC_ICGC
#     tissues_to_consider = config.tissues_to_consider
#     tss_fragment_filter = config.tss_fragment_filter
#     ML_model = config.ML_model
#     lung_subtyped = config.lung_subtyped
#     woo_pcawg = config.woo_pcawg
#     histologically_subtyped_mutations = config.histologically_subtyped_mutations
#     de_novo_seurat_clustering = config.de_novo_seurat_clustering
#     meso = config.meso
#     RNA_subtyped = config.RNA_subtyped
#     per_donor = config.per_donor
#     donor_range = config.donor_range
#     test_backward_selection_iters = config.test_backward_selection_iters
#     seed_range = config.seed_range
#     start, end = map(int, seed_range.split('-'))
#     seed_range = range(start, end + 1)
#     n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
#     n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
#     top_features_to_plot = config.top_features_to_plot
#     save_test_set_perf = config.save_test_set_perf
#     make_plots = config.make_plots
#
#     return {
#         'cancer_types': cancer_types,
#         'datasets': datasets,
#         'scATAC_cell_number_filter': scATAC_cell_number_filter,
#         'annotation_dir': annotation_dir,
#         'SCLC': SCLC,
#         'CPTAC': CPTAC,
#         'combined_CPTAC_ICGC': combined_CPTAC_ICGC,
#         'tissues_to_consider': tissues_to_consider,
#         'tss_fragment_filter': tss_fragment_filter,
#         'ML_model': ML_model,
#         'lung_subtyped': lung_subtyped,
#         'woo_pcawg': woo_pcawg,
#         'histologically_subtyped_mutations': histologically_subtyped_mutations,
#         'de_novo_seurat_clustering': de_novo_seurat_clustering,
#         'meso': meso,
#         'RNA_subtyped': RNA_subtyped,
#         'per_donor': per_donor,
#         'donor_range': donor_range,
#         'test_backward_selection_iters': test_backward_selection_iters,
#         'seed_range': seed_range,
#         'n_optuna_trials_prebackward_selection': n_optuna_trials_prebackward_selection,
#         'n_optuna_trials_backward_selection': n_optuna_trials_backward_selection,
#         'top_features_to_plot': top_features_to_plot,
#         'save_test_set_perf': save_test_set_perf,
#         'make_plots': make_plots
#     }

# cancer_types = config.cancer_types
# # run_clustered_mutations = config.clustered_mutations
# datasets = sorted(config.datasets)
# scATAC_cell_number_filter = config.scATAC_cell_number_filter
# annotation_dir = config.annotation_dir
# SCLC = config.SCLC
# CPTAC = config.CPTAC
# combined_CPTAC_ICGC = config.combined_CPTAC_ICGC
# tissues_to_consider = config.tissues_to_consider
# # tss_filtered = config.tss_filtered
# tss_fragment_filter = config.tss_fragment_filter
# ML_model = config.ML_model
# lung_subtyped = config.lung_subtyped
# woo_pcawg = config.woo_pcawg
# histologically_subtyped_mutations = config.histologically_subtyped_mutations
# de_novo_seurat_clustering = config.de_novo_seurat_clustering
# meso = config.meso
# RNA_subtyped = config.RNA_subtyped
# per_donor = config.per_donor
# donor_range = config.donor_range
# test_backward_selection_iters = config.test_backward_selection_iters
# seed_range = config.seed_range
# start, end = map(int, seed_range.split('-'))
# seed_range = range(start, end + 1)
# n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
# n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
# top_features_to_plot = config.top_features_to_plot
# save_test_set_perf = config.save_test_set_perf
# make_plots = config.make_plots
# iters_dont_skip = config.iters_dont_skip
# num_iter_skips = config.num_iter_skips