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
    parser.add_argument('--datasets', nargs="+", type=str,
                        help='which sc-ATACseq datasets to analyze', required=True)
    parser.add_argument('--scATAC_cell_number_filter', type=int,
                        help='minimum number of cells per cell type in scATAC', default=100)
    parser.add_argument('--annotation_dir', type=str,
                        help='name of annotation directory', default="default_annotation")
    parser.add_argument('--tss_fragment_filter', nargs="+", type=str,
                        help='tss fragment filters to consider', default=None)
    parser.add_argument('--tissues_to_consider', nargs="+", type=str, default="all")
    parser.add_argument("--ML_model", type=str, default="XGB")
    parser.add_argument('--test_backward_selection_iters', type=int, nargs="+", default=None)
    parser.add_argument('--seed_range', type=str)
    parser.add_argument('--fold_for_test_set', type=int)
    parser.add_argument('--n_optuna_trials_prebackward_selection', type=int, default=None)
    parser.add_argument('--n_optuna_trials_backward_selection', type=int, default=None)
    parser.add_argument('--top_features_to_plot', nargs="+", type=int)
    parser.add_argument("--save_test_set_perf", action="store_true", default=False)
    parser.add_argument("--make_plots", action="store_true", default=False)
    parser.add_argument("--feature_importance_method", type=str, default="permutation_importance")
    parser.add_argument("--sqlite", action="store_true", default=False)
    parser.add_argument('--test_set_perf_num_features', nargs="+", type=int)
    parser.add_argument('--debug_bfs', action="store_true", default=False)
    parser.add_argument('--donor_range', type=range_type, help='Specify a range in the format start-end',
                        default=None)
    parser.add_argument("--per_donor", action="store_true", default=False)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--SCLC", action="store_true", default=False)
    group.add_argument("--lung_subtyped", action="store_true", default=False)
    group.add_argument("--woo_pcawg", action="store_true", default=False)
    group.add_argument("--histologically_subtyped_mutations", action="store_true", default=False)
    group.add_argument("--de_novo_seurat_clustering", action="store_true", default=False)
    group.add_argument("--CPTAC", action="store_true", default=False)
    group.add_argument("--meso", action="store_true", default=False)
    group.add_argument("--hundred_kb", action="store_true", default=False)
    group.add_argument("--expanded_hundred_kb", action="store_true", default=False)
    group.add_argument("--combined_CPTAC_ICGC", action="store_true", default=False)
    group.add_argument("--RNA_subtyped", action="store_true", default=False)
    return parser
