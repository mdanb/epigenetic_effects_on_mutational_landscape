from ML_utils import *
import numpy as np
import matplotlib.pyplot as plt
from config import create_parser

parser = create_parser()
parser.add_argument('--error_analysis', action="store_true", default=False)
parser.add_argument('--count_bin_sums', action="store_true", default=False)

config = parser.parse_args()

cancer_types = config.cancer_types
datasets = sorted(config.datasets)
scATAC_cell_number_filter = config.scATAC_cell_number_filter
annotation_dir = config.annotation_dir
SCLC = config.SCLC
CPTAC = config.CPTAC
combined_CPTAC_ICGC = config.combined_CPTAC_ICGC
ML_model = config.ML_model
lung_subtyped = config.lung_subtyped
woo_pcawg = config.woo_pcawg
histologically_subtyped_mutations = config.histologically_subtyped_mutations
de_novo_seurat_clustering = config.de_novo_seurat_clustering
meso = config.meso
RNA_subtyped = config.RNA_subtyped
per_donor = config.per_donor
donor_range = config.donor_range
seed_range = config.seed_range
start, end = map(int, seed_range.split('-'))
seed_range = range(start, end + 1)
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
feature_importance_method = config.feature_importance_method
tss_filter = config.tss_fragment_filter
num_features = config.error_analysis_num_features
fold_for_test_set = config.fold_for_test_set
error_analysis = config.error_analysis
count_bin_sums = config.count_bin_sums

if error_analysis:
    for seed in seed_range:
        for cancer_type in cancer_types:
            scATAC_df, cancer_specific_mutations = load_data(meso, SCLC, lung_subtyped, woo_pcawg,
                                                              histologically_subtyped_mutations,
                                                              de_novo_seurat_clustering,
                                                              CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
                                                              datasets, scATAC_cell_number_filter,
                                                              annotation_dir, cancer_type)
            scATAC_sources = construct_scATAC_sources(datasets)
            scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir, seed,
                                              fold_for_test_set)

            model = load_n_features_backwards_elim_models(num_features, scATAC_df.shape[1], cancer_type, ML_model,
                                                          scATAC_dir, feature_importance_method, full_data_trained=True)
            top_features = model.feature_names_in_
            scATAC_df = scATAC_df.loc[:, top_features]
            # X_train, _, y_train, _ = get_train_test_split(scATAC_df, cancer_specific_mutations, 0.10, seed)
            errors = apply_func_to_kfolds(scATAC_df, cancer_specific_mutations, compute_error, estimator=model,
                                          train_new_model=True)
            abs_errors = []
            percent_errors = []
            preds = []
            for fold_results in errors:
                abs_errors.append(fold_results["abs_err"])
                percent_errors.append(fold_results["percent_err"])
                preds.append(fold_results["preds"])
            abs_errors = pd.concat(abs_errors)
            percent_errors = pd.concat(percent_errors)
            preds = np.concatenate(preds)
            # errors = [err for err_list in errors for err in err_list]
            fig, ax = plt.subplots(5, 1, figsize=(15,15), sharex=True)
            x = np.arange(len(abs_errors))
            to_plot = [abs_errors, percent_errors, cancer_specific_mutations, preds, scATAC_df.to_numpy().reshape(-1)]
            ylabel_list = ["Absolute error", "Percent error", "Mutation counts", "Predictions", "scATAC"]
            for i in range(5):
                ax[i].bar(x, np.asarray(to_plot[i], dtype=np.float32), width=4)
                # for j in range(0, 2200, 213):
                #     ax[i].axvline(x=j, c="red")
                ax[i].set_ylabel(ylabel_list[i])
            # plt.title("Percent error versus chromosomal bins")
            ax[-1].set_xticks(np.arange(min(x), max(x)+1, 100))
            ax[-1].set_xlabel("Bin")
            for a in ax:
                a.tick_params(axis='x', labelsize=8)
            plt.savefig(f"../../figures/models/{ML_model}/{cancer_type}/{scATAC_dir}/errors.png")

# Use only seed=1 since splits are the same for all seeds
if count_bin_sums:
    for cancer_type in cancer_types:
        for i in range(1, 10):
            scATAC_sources = construct_scATAC_sources(datasets)
            scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir,
                                              seed=1, fold_for_test_set=i)
            top_feature = open(f"models/{ML_model}/{cancer_type}/{scATAC_dir}/backwards_elimination_results/"
                               f"top_features_iteration_19_by_permutation_importance.txt", "r").readline().strip()[3:]

            X_test = pd.read_csv(f"models/{ML_model}/{cancer_type}/{scATAC_dir}/X_test.csv", index_col=0)
            top_feature_x_test = X_test.loc[:, top_feature]
            top_feat_total_x_test_scatac_counts = sum(top_feature_x_test)
            X_train = pd.read_csv(f"models/{ML_model}/{cancer_type}/{scATAC_dir}/X_train.csv", index_col=0)
            top_feature_x_train = X_train.loc[:, top_feature]
            top_feat_total_x_train_scatac_counts = sum(top_feature_x_train)
            y_test = pd.read_csv(f"models/{ML_model}/{cancer_type}/{scATAC_dir}/y_test.csv", index_col=0)
            total_y_test_counts = y_test.sum()[0]
            y_train = pd.read_csv(f"models/{ML_model}/{cancer_type}/{scATAC_dir}/y_train.csv", index_col=0)
            total_y_train_counts = y_train.sum()[0]

