from ML_utils import *
import numpy as np
# import matplotlib.pyplot as plt
from config import create_parser
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
import os.path as osp
# import pyreadr

parser = create_parser()
parser.add_argument('--bins_error_analysis', action="store_true", default=False)
parser.add_argument('--count_bin_sums', action="store_true", default=False)
parser.add_argument('--error_analysis_num_features', type=int, default=1)

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
hundred_kb = config.hundred_kb
expanded_hundred_kb = config.expanded_hundred_kb
per_donor = config.per_donor
donor_range = config.donor_range
# seed_range = config.seed_range
# start, end = map(int, seed_range.split('-'))
# seed_range = range(start, end + 1)
n_optuna_trials_prebackward_selection = config.n_optuna_trials_prebackward_selection
n_optuna_trials_backward_selection = config.n_optuna_trials_backward_selection
feature_importance_method = config.feature_importance_method
tss_filter = config.tss_fragment_filter
error_analysis_num_features = config.error_analysis_num_features
fold_for_test_set = config.fold_for_test_set
bins_error_analysis = config.bins_error_analysis
count_bin_sums = config.count_bin_sums
aggregated_per_donor = config.aggregated_per_donor
tissues_to_consider = config.tissues_to_consider
grid_analysis = config.grid_analysis
grid_cell_types = config.grid_cell_types
cell_types_keep = config.cell_types_keep
hierarchically_subtyped_mutations = config.hierarchically_subtyped_mutations
mm = config.mm
msi_high = config.msi_high

subsampled_mutations = config.subsampled_mutations

def conduct_test(errors, underestimated=None, two_sided=None):
    assert not ((underestimated is not None and two_sided is not None) or
                (underestimated is None and two_sided is None))
    normalized_errors = zscore(errors)
    if two_sided:
        p_values = 2 * (1 - norm.cdf(abs(normalized_errors)))
    else:
        if underestimated:
            p_values = 1 - norm.cdf(normalized_errors)
        else:
            p_values = norm.cdf(normalized_errors)
    rejected, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    return rejected, q_values, normalized_errors, p_values

if bins_error_analysis:
#     for seed in seed_range:
    scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir,
                                    hundred_kb, expanded_hundred_kb, tissues_to_consider, grid_analysis,
                                    grid_cell_types, cell_types_keep)
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]

    for cancer_type in cancer_types:
        mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                      histologically_subtyped_mutations, de_novo_seurat_clustering,
                                      CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, cancer_type,
                                      hundred_kb, expanded_hundred_kb, aggregated_per_donor,
                                      hierarchically_subtyped_mutations, mm, msi_high, subsampled_mutations)
        scATAC_df, cancer_specific_mutations = prep_and_align_mutations_with_scatac(scATAC_df,
                                                                                    mutations_df,
                                                                                    cancer_type,
                                                                                    hundred_kb,
                                                                                    expanded_hundred_kb,
                                                                                    per_donor,
                                                                                    subsampled_mutations,
                                                                                    None)

        # scATAC_df, cancer_specific_mutations = load_data(meso, SCLC, lung_subtyped, woo_pcawg,
        #                                                   histologically_subtyped_mutations,
        #                                                   de_novo_seurat_clustering,
        #                                                   CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
        #                                                   datasets, scATAC_cell_number_filter,
        #                                                   annotation_dir, cancer_type, hundred_kb)
        scATAC_sources = construct_scATAC_sources(datasets)
        multiseed_errors = []
        # multiseed_percent_errors = []
        multiseed_preds = []
        actual = []
        atac = []
        for seed in range(1, 11):
            abs_errors = []
            # percent_errors = []
            preds = []
            for fold_for_test_set in range(0, 10):
                scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir,
                                                  hundred_kb, expanded_hundred_kb, tissues_to_consider,
                                                  seed=seed, fold_for_test_set=fold_for_test_set,
                                                  all_seeds=False)
                model = load_n_features_backwards_elim_models(error_analysis_num_features, scATAC_df.shape[1],
                                                              cancer_type, ML_model, scATAC_dir, feature_importance_method,
                                                              full_data_trained=True)
                top_features = model.feature_names_in_
                _, X_test, _, y_test = get_train_test_split(scATAC_df, cancer_specific_mutations, 10, fold_for_test_set)
                x = X_test.loc[:, top_features]
                pred = model.predict(x)
                error = y_test - pred
                # errors = compute_error(X_test.loc[:, top_features], y_test, estimator=model, train_new_model=False,
                #                        i=fold_for_test_set)
                # errors = apply_func_to_kfolds(scATAC_df.loc[:, top_features], cancer_specific_mutations, compute_error,
                #                               estimator=model,
                #                               train_new_model=False)
                # for fold_results in errors:
                if seed == 1:
                    actual.append(y_test)
                    atac.append(x.squeeze())
                abs_errors.append(error)
                # percent_errors.append(errors["percent_err"])
                preds.append(pred)
            multiseed_errors.append(pd.concat(abs_errors))
            # multiseed_percent_errors.append(pd.concat(percent_errors))
            multiseed_preds.append(np.concatenate(preds))

        multiseed_errors = pd.concat(multiseed_errors, axis=1)
        avg_multiseed_errors = multiseed_errors.mean(axis=1)
        avg_errors_se = multiseed_errors.sem(axis=1)
        median_multiseed_errors = multiseed_errors.median(axis=1)
        # multiseed_percent_errors = pd.concat(multiseed_percent_errors, axis=1).mean(axis=1)
        avg_multiseed_preds = np.mean(multiseed_preds, axis=0)
        median_multiseed_preds = np.median(multiseed_preds, axis=0)
        actual = np.concatenate(actual)
        atac = np.concatenate(atac)
        # multiseed_percent_errors.plot.hist()
        # percent_errors_rejected, percent_errors_q_values, \
        # normalized_percent_errors, percent_errors_p_values = conduct_test(multiseed_percent_errors, two_sided=True)
        # absolute_errors_rejected_under, absolute_errors_q_values_under, \
        # normalized_absolute_errors, absolute_errors_p_values_under = conduct_test(multiseed_errors,
        #                                                                                 underestimated=True)
        # absolute_errors_rejected_over, absolute_errors_q_values_over, \
        # _, absolute_errors_p_values_over = conduct_test(multiseed_errors, underestimated=True)

        # normalized_percent_errors = zscore(multiseed_percent_errors)
        # p_values = 2 * (1 - norm.cdf(abs(normalized_percent_errors)))
        # rejected, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

        df = pd.DataFrame({
                           # "percent_errors_p-value": percent_errors_p_values,
                           # "normalized_percent_error": normalized_percent_errors,
                           # "percent_error_q-value": percent_errors_q_values,
                           # "percent_error_rejected": percent_errors_rejected,
                           "avg_absolute_error": avg_multiseed_errors,
                           "avg_errors_se": avg_errors_se,
                           "median_absolute_error": median_multiseed_errors,
                           "avg_prediction": avg_multiseed_preds,
                           "median_prediction": median_multiseed_preds,
                           "actual": actual,
                           "atac": atac
                           # "avg_percent_error": multiseed_percent_errors,
                           # "normalized_absolute_error": normalized_absolute_errors,
                           # "absolute_errors_p-value_under": absolute_errors_p_values_under,
                           # "absolute_error_q-value_under": absolute_errors_q_values_under,
                           # "absolute_error_rejected_under": absolute_errors_rejected_under,
                           # "absolute_error_rejected_over": absolute_errors_rejected_over,
                           # "absolute_errors_p-value_over": absolute_errors_p_values_over,
                           # "absolute_error_q-value_over": absolute_errors_q_values_over,
                           })

        fp = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir,
                                  hundred_kb, expanded_hundred_kb, tissues_to_consider, all_seeds=True)
        # Good style / Cool pattern
        current_fp = osp.dirname(osp.realpath(__file__))
        fp = osp.join(current_fp, "..", "..", "figures", "models", ML_model,
                      cancer_type, fp)
        multiseed_errors.to_csv(osp.join(fp, "errors_df.csv"))

        # df = df.sort_values(by="absolute_error_q-value")
        # if hundred_kb:
        #     ranges = pyreadr.read_r(osp.join(current_fp, "..", "..", "data", "100kb_interval_ranges.Rdata"))
        # else:
        #     ranges = pyreadr.read_r(osp.join(current_fp, "..", "..", "data", "mutation_data",
        #                                      "hg19.1Mb.ranges.Polak.Nature2015.RData"))
        df.to_csv(osp.join(fp, "errors_df_summary.csv"))

        # multiseed_preds = np.average(multiseed_preds, axis=0)


            # abs_errors = pd.concat(abs_errors)
            # percent_errors = pd.concat(percent_errors)
            # preds = np.concatenate(preds)
        # # errors = [err for err_list in errors for err in err_list]
        # fig, ax = plt.subplots(5, 1, figsize=(15,15), sharex=True)
        # x = np.arange(len(abs_errors))
        # to_plot = [abs_errors, percent_errors, cancer_specific_mutations, preds, scATAC_df.to_numpy().reshape(-1)]
        # ylabel_list = ["Absolute error", "Percent error", "Mutation counts", "Predictions", "scATAC"]
        # for i in range(5):
        #     ax[i].bar(x, np.asarray(to_plot[i], dtype=np.float32), width=4)
        #     # for j in range(0, 2200, 213):
        #     #     ax[i].axvline(x=j, c="red")
        #     ax[i].set_ylabel(ylabel_list[i])
        # # plt.title("Percent error versus chromosomal bins")
        # ax[-1].set_xticks(np.arange(min(x), max(x)+1, 100))
        # ax[-1].set_xlabel("Bin")
        # for a in ax:
        #     a.tick_params(axis='x', labelsize=8)
        # plt.savefig(f"../../figures/models/{ML_model}/{cancer_type}/{scATAC_dir}/errors.png")

# Use only seed=1 since splits are the same for all seeds
if count_bin_sums:
    for cancer_type in cancer_types:
        for i in range(1, 11):
            scATAC_sources = construct_scATAC_sources(datasets)
            scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir,
                                               hundred_kb, expanded_hundred_kb, seed=1, fold_for_test_set=i)
            top_feature = open(osp.join(ML_model, cancer_type,scATAC_dir, 'backwards_elimination_results',
                               'top_features_iteration_19_by_permutation_importance.txt'), "r").readline().strip()[3:]
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

