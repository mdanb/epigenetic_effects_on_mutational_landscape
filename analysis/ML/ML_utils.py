import pyreadr
from natsort import natsorted
from itertools import chain
import pandas as pd
import numpy as np
import pickle
import os
from sklearn.metrics import r2_score
from xgboost import XGBRegressor
from sklearn.model_selection import cross_validate, KFold
import optuna
import subprocess
from sklearn.inspection import permutation_importance
from pathlib import Path
import re

### Load Data helpers ###
def load_data(meso, SCLC, lung_subtyped, woo_pcawg,
              histologically_subtyped_mutations, de_novo_seurat_clustering,
              CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor,
              datasets, scATAC_cell_number_filter, annotation_dir,
              cancer_type_or_donor_id, tss_filter=None):
    mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                  histologically_subtyped_mutations, de_novo_seurat_clustering,
                                  CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, cancer_type_or_donor_id)

    scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir)
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    chr_keep = pd.read_csv("../../data/processed_data/chr_keep.csv", index_col=0)
    mutations_df = mutations_df.loc[chr_keep["chr"]]
    if not pd.isna(mutations_df).any().any():
        # for compatibility when chr ranges only include chr_keep but not others
        mutations_df = add_na_ranges(mutations_df)
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)
    cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type_or_donor_id)
    return scATAC_df, cancer_specific_mutations


def load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                   histologically_subtyped_mutations, de_novo_seurat_clustering,
                   CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, cancer_type):
    if meso:
        mutations_df = load_meso()
    elif SCLC:
        mutations_df = load_sclc_mutations()
    elif lung_subtyped:
        mutations_df = load_subtyped_lung_mutations()
    elif woo_pcawg:
        mutations_df = load_woo_pcawg_mutations()
    elif histologically_subtyped_mutations:
        mutations_df = load_histologically_subtyped_mutations()
    elif de_novo_seurat_clustering:
        mutations_df = load_de_novo_seurat_clustered_cancers(cancer_type)
    elif CPTAC:
        mutations_df = load_CPTAC()
    elif combined_CPTAC_ICGC:
        mutations_df = load_combined_CPTAC_ICGC()
    elif RNA_subtyped:
        mutations_df = load_RNA_subtyped_mutations()
    elif per_donor:
        mutations_df = load_per_donor_mutations()
    else:
        mutations_df = load_agg_mutations()
    return mutations_df


def load_scATAC(scATAC_path):
    scATAC_df = pyreadr.read_r(scATAC_path)
    scATAC_df = scATAC_df[None]
    scATAC_df = scATAC_df.T
    return scATAC_df


def load_scATAC_metadata(metadata_path):
    metadata = pyreadr.read_r(metadata_path)
    metadata = metadata[None]
    return metadata


def load_per_donor_mutations(cancer_type):
    df = pd.read_csv(f"../../data/processed_data/per_patient_mutations/{cancer_type}_per_donor.csv",
                     index_col=0)
    return df.loc[natsorted(df.index)]


def load_CPTAC():
    df = pd.read_csv(f"../../data/processed_data/CPTAC_mutations.csv",
                     index_col=0)
    return df.loc[natsorted(df.index)]


def load_combined_CPTAC_ICGC():
    df = pd.read_csv(f"../../data/processed_data/combined_CPTAC_ICGC_mutations.csv",
                     index_col=0)
    return df.loc[natsorted(df.index)]


def load_agg_mutations():
    df = pd.read_csv("../../data/processed_data/mut_count_data.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_woo_pcawg_mutations():
    df = pd.read_csv("../../data/processed_data/pcawg_agg_woo.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_RNA_subtyped_mutations():
    df = pd.read_csv("../../data/processed_data/RNA_subtyped_cancers.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_meso():
    df = pd.read_csv("../../data/processed_data/mesothelioma.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_sclc_mutations():
    df = pd.read_csv("../../data/processed_data/sclc_count_overlaps.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_subtyped_lung_mutations():
    df = pd.read_csv("../../data/processed_data/lung_subtyped_mutation.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_histologically_subtyped_mutations():
    df = pd.read_csv("../../data/processed_data/histologically_subtyped_mutations.csv",
                       index_col=0)
    return df.loc[natsorted(df.index)]


def load_de_novo_seurat_clustered_cancers(cancer_type):
    print("Loading De Novo Seurat clustered cancers")
    seurat_cluster_settings = cancer_type.split("x")[1]
    cancer_type = cancer_type.split("x")[0]
    df = pd.read_csv(f"../../data/processed_data/de_novo_seurat_clustered_mutations/{cancer_type}/"
                     f"{seurat_cluster_settings}.csv", index_col=0)
    return df.loc[natsorted(df.index)]


#### Filter Data Helpers ####
def filter_agg_data(scATAC_df, mutations_df):
    idx_select = ~pd.isnull(mutations_df).any(axis=1)
    scATAC_df = scATAC_df.loc[idx_select]
    mutations_df = mutations_df[idx_select]
    return scATAC_df, mutations_df


def filter_mutations_by_cancer(mutations_df, cancer_type):
    cancer_type_specific_mut_idx = [i for i, s in enumerate(mutations_df.columns.values) if cancer_type == s][0]
    cancer_specific_mutations = mutations_df.iloc[:, cancer_type_specific_mut_idx]
    return cancer_specific_mutations


def filter_clustered_data(scATAC_df, mutations_df):
    return pd.DataFrame(scATAC_df, index=mutations_df.index)


def filter_scATAC_df_by_num_cell_per_cell_type(scATAC_df, scATAC_cell_number_filter, metadata):
    metadata = metadata.loc[metadata["num_cells"].astype(int) >= scATAC_cell_number_filter, :]
    try:
        keep = [tissue + " " + cell_type for tissue, cell_type in zip(metadata["tissue_name"], metadata["cell_type"])]
    except KeyError:
        keep = [tissue + " " + cell_type for tissue, cell_type in zip(metadata["tissue"], metadata["cell_type"])]
    scATAC_df = scATAC_df.loc[:, scATAC_df.columns.isin(keep)]
    return scATAC_df

#### Dataframe curation helpers ####
def add_na_ranges(mutations_df):
    full_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
    mutations_df = pd.merge(full_ranges, mutations_df, left_on = "x", right_index=True,
                            how="outer").set_index("x")
    return mutations_df.loc[natsorted(mutations_df.index)]


def add_dataset_origin_to_cell_types(df, dataset):
    if dataset == "Bingren":
        df.columns = [c + " BR" for c in df.columns]
    elif dataset == "Shendure":
        df.columns = [c + " SH" for c in df.columns]
    elif dataset == "Tsankov":
        df.columns = [c + " TS" for c in df.columns]
    elif dataset == "Greenleaf_brain":
        df.columns = [c + " GL_Br" for c in df.columns]
    elif dataset == "Greenleaf_pbmc_bm":
        df.columns = [c + " GL_BlBm" for c in df.columns]
    elif dataset == "Yang_kidney":
        df.columns = [c + " Y_K" for c in df.columns]
    elif dataset == "Rawlins_fetal_lung":
        df.columns = [c + " R_Fl" for c in df.columns]
    elif dataset == "Greenleaf_colon":
        df.columns = [c + " GL_Co" for c in df.columns]
    elif dataset == "Wang_lung":
        df.columns = [c + " W_L" for c in df.columns]
    return df


def construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir):
    datasets_combined_count_overlaps = []
    for dataset in datasets:
        if tss_filter:
            print(f"Loading TSS filtered scATAC from {dataset}...")
            tss_filtered_root = "../../data/processed_data/count_overlap_data/tsse_filtered"
            chr_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
            scATAC_df = load_scATAC(f"{tss_filtered_root}/{dataset}/combined/{annotation_dir}/" \
                                    f"combined_{tss_filter}_fragments.rds").T
            print("Loaded!")
            scATAC_df.index = chr_ranges["x"].values
            datasets_combined_count_overlaps.append(scATAC_df)
        else:
            print(f"Loading scATAC from {dataset}...")
            scATAC_df = load_scATAC("../../data/processed_data/count_overlap_data/combined_count_overlaps" \
            f"/{annotation_dir}/{dataset}_combined_count_overlaps.rds")
            print("Loaded!")
            metadata = load_scATAC_metadata("../../data/processed_data/count_overlap_data/combined_count_overlaps" \
            f"/{annotation_dir}/{dataset}_combined_count_overlaps_metadata.rds")
            scATAC_df = filter_scATAC_df_by_num_cell_per_cell_type(scATAC_df, scATAC_cell_number_filter, metadata)
            datasets_combined_count_overlaps.append(scATAC_df)

    for idx, dataset in enumerate(datasets):
        datasets_combined_count_overlaps[idx] = [add_dataset_origin_to_cell_types(datasets_combined_count_overlaps[idx],
                                                 dataset)]
    scATAC_df = pd.concat(chain(*datasets_combined_count_overlaps), axis=1)
    return scATAC_df


#### Split Train/Test helpers ####
def get_train_test_split(X, y, n_splits, fold_for_test_set):
    # Don't shuffle, to make problem harder as explained elsewhere
    kf = KFold(n_splits=n_splits, shuffle=False)
    train_index = -1
    test_index = -1
    for i, (train_idx, test_idx) in enumerate(kf.split(X)):
        if fold_for_test_set == i:
            train_index = train_idx
            test_index = test_idx
            break
    X_train, y_train = X.iloc[train_index], y.iloc[train_index]
    X_test, y_test = X.iloc[test_index], y.iloc[test_index]
    return X_train, X_test, y_train, y_test


#### Feature Selection helpers ####
def calculate_permutation_importance(X_valid, y_valid, **kwargs):
    estimator = kwargs["estimator"]
    n_repeats = kwargs["n_repeats"]
    seed = kwargs["seed"]
    if kwargs["i"] is not None:
        estimator = estimator[kwargs["i"]]
        print(f"Calculating permutation importance for fold {kwargs['i']}...")
    out = permutation_importance(estimator, X_valid, y_valid, n_repeats=n_repeats, random_state=seed, n_jobs=-1)
    return out


def get_top_n_features(best_model_fulldatatrained, best_model_perfoldtrained: list,
                       n, features, feature_importance_method, X, y, seed,
                       df_save=None, fp_for_fi=None, best_cv_score=None):
    if df_save is not None and n + 1 in df_save["num_features"].array:
        feature_importances = df_save.loc[df_save["num_features"] == n + 1][feature_importance_method]
    else:
        # kf = KFold(n_splits=10)
        std = [np.NaN] * len(features)
        if feature_importance_method == "default_importance":
            feature_importances = best_model_fulldatatrained.feature_importances_
        elif feature_importance_method == "permutation_importance":
            assert not (X is None or y is None or seed is None)
            print("Computing permutation importance per fold...")
            func_out = apply_func_to_kfolds(X, y, calculate_permutation_importance, n_repeats=10,
                                                     seed=seed, estimator=best_model_perfoldtrained)
            feature_importance_means = []
            std = []
            for fold_results in func_out:
                feature_importance_means.append(fold_results["importances_mean"])
                std.append(fold_results["importances_std"])
            feature_importances = np.stack(feature_importance_means).mean(axis=0)
            std = np.stack(std).mean(axis=0)
            print("Done computing feature importances!")
        if df_save is not None:
            df_curr = pd.DataFrame((features, feature_importances, [len(features)] * len(features),
                                    [best_cv_score] * len(features), std)).T
            df_curr.columns = df_save.columns
            df_save = pd.concat((df_save, df_curr))
            df_save.to_csv(fp_for_fi, index=False)
    feat_importance_idx = np.argsort(feature_importances)[::-1]
    top_n_feats = features[feat_importance_idx][:n]
    return top_n_feats, df_save


def print_and_save_features(features, filepath, top=True):
    n = len(features)
    if top:
        print(f"Top {n} features")
    f = open(filepath, "w")
    for idx, feature in enumerate(features):
        if top:
            print(f"{idx+1}. {feature}")
            f.write(f"{idx+1}. {feature}\n")
        else:
            print(f"-{feature}")
            f.write(f"-{feature}\n")


def backward_eliminate_features(X_train, y_train, backwards_elim_dir,
                                ML_model, scATAC_dir, cancer_type_or_donor_id, seed,
                                n_optuna_trials_backward_selection, feature_importance_method, sqlite,
                                best_model_fulldatatrained, best_model_perfoldtrained,
                                starting_n):
    figure_path = os.path.join("../../figures/models",
                               ML_model,
                               cancer_type_or_donor_id,
                               scATAC_dir,
                               "backwards_elimination_results")
    Path(figure_path).mkdir(parents=True, exist_ok=True)
    filename_df_for_fi = "df_for_feature_importance_plots"
    if feature_importance_method != "default_importance":
        filename_df_for_fi = filename_df_for_fi + f"_{feature_importance_method}"
    filename_df_for_fi = filename_df_for_fi + ".csv"
    fp_for_fi = figure_path + "/" + filename_df_for_fi
    print(f"starting_n is {starting_n}")
    if starting_n is not None:
        if not os.path.exists(f"{backwards_elim_dir}/all_features_rankings_by_{feature_importance_method}.txt"):
            ranked_features, df_save = get_top_n_features(best_model_fulldatatrained,
                                                           best_model_perfoldtrained,
                                                           len(X_train.columns.values),
                                                           X_train.columns.values,
                                                           feature_importance_method,
                                                           X_train, y_train, seed)
            print_and_save_features(ranked_features, filepath=f"{backwards_elim_dir}/"
                                    f"all_features_rankings_by_{feature_importance_method}.txt",
                                    top=True)
        else:
            with open(f"{backwards_elim_dir}/all_features_rankings_by_{feature_importance_method}.txt", "r") as f:
                ranked_features = f.read().splitlines()
                ranked_features = [re.sub("^\d+\.\s*", "", feature) for feature in ranked_features]
        top_n_features = ranked_features[:starting_n]
        X_train = X_train.loc[:, top_n_features]
        print_and_save_features(top_n_features,
                                filepath=f"{backwards_elim_dir}/starter_model_top_features_by_"
                                         f"{feature_importance_method}.txt",
                                top=True)
        num_iterations = len(top_n_features)

    else:
        print_and_save_features(X_train.columns,
                                filepath=f"{backwards_elim_dir}/starter_features.txt",
                                top=False)
        num_iterations = X_train.shape[1]

    if os.path.exists(fp_for_fi):
        df_save = pd.read_csv(fp_for_fi)
    else:
        df_save = pd.DataFrame(columns=["features", feature_importance_method, "num_features", "score", "std"])

    for idx in range(1, num_iterations + 1):
        print(f"Running BFS iteration {idx}...")
        model_optimizer = ModelOptimizer(backwards_elim_dir + "/" + f"model_optimizer_iteration_{idx}.pkl")
        top_features_filepath = f"{backwards_elim_dir}/top_features_iteration_{idx}"
        model_savefile = f"model_iteration_{idx}"
        per_fold_model_savefile = f"per_fold_model_iteration_{idx}"

        study_name = f"{cancer_type_or_donor_id}_iter_{idx}_{scATAC_dir}"
        if feature_importance_method != "default_importance":
            study_name = study_name + f"_feature_importance_{feature_importance_method}"
            top_features_filepath = top_features_filepath + f"_by_{feature_importance_method}"
            model_savefile = model_savefile + f"_feature_importance_{feature_importance_method}"
            per_fold_model_savefile = per_fold_model_savefile + f"_feature_importance_{feature_importance_method}"

        top_features_filepath = top_features_filepath + ".txt"
        model_savefile = model_savefile + ".pkl"
        per_fold_model_savefile = per_fold_model_savefile + ".pkl"

        study = model_optimizer.optimize_optuna_study(study_name=study_name, ML_model=ML_model, X_train=X_train,
                                                      y_train=y_train, seed=seed,
                                                      n_optuna_trials=n_optuna_trials_backward_selection,
                                                      sqlite=sqlite)
        best_model_perfoldtrained = model_optimizer.best_model_perfoldtrained
        if not os.path.exists(f"{backwards_elim_dir}/{model_savefile}"):
            best_params = study.best_params
            if ML_model == "XGB":
                best_model_fulldatatrained = XGBRegressor(**best_params)

            best_model_fulldatatrained.fit(X=X_train, y=y_train)
            pickle.dump(best_model_fulldatatrained, open(f"{backwards_elim_dir}/{model_savefile}", 'wb'))
        else:
            best_model_fulldatatrained = pickle.load(open(f"{backwards_elim_dir}/{model_savefile}", "rb"))

        if not os.path.exists(f"{backwards_elim_dir}/{per_fold_model_savefile}"):
            pickle.dump(best_model_perfoldtrained, open(f"{backwards_elim_dir}/{per_fold_model_savefile}", 'wb'))
        else:
            best_model_perfoldtrained = pickle.load(open(f"{backwards_elim_dir}/{per_fold_model_savefile}", "rb"))

        best_trial = study.best_trial
        best_cv_score = best_trial.value

        # # No feature elimination if only 1 feature is left
        # if idx != num_iterations:
        top_n_feats, df_save = get_top_n_features(best_model_fulldatatrained,
                                                   best_model_perfoldtrained,
                                                   len(X_train.columns.values) - 1,
                                                   X_train.columns.values,
                                                   feature_importance_method,
                                                   X_train, y_train, seed,
                                                   df_save=df_save, fp_for_fi=fp_for_fi,
                                                   best_cv_score=best_cv_score)
        removed_feature = X_train.columns[~(X_train.columns.isin(top_n_feats))].item()
        print(f"Removing: {removed_feature}")
        X_train = X_train.loc[:, top_n_feats]
        # Don't save 0 features
        if idx != num_iterations:
            print_and_save_features(top_n_feats, filepath=top_features_filepath, top=True)


#### Model train/val/test helpers ####
class ModelOptimizer:
    def __init__(self, savepath):
        self.best_model_perfoldtrained = None
        self._best_score = float('-inf')
        self._savepath = savepath
        if os.path.exists(savepath):
            model_optimizer = pickle.load(open(savepath, "rb"))
            self._best_score = model_optimizer["best_score"]
            self.best_model_perfoldtrained = model_optimizer["best_model_perfoldtrained"]

    def _optuna_objective(self, trial, ML_model, X, y, seed):
        if ML_model == "XGB":
            param = {
                'n_estimators': trial.suggest_int('n_estimators', 100, 500),
                'max_depth': trial.suggest_int('max_depth', 3, 10),
                'learning_rate': trial.suggest_float('learning_rate', 1e-8, 1.0, log=True),
                'subsample': trial.suggest_float('subsample', 0.1, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 6),
                'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
                'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
            }
            model = XGBRegressor(**param, random_state=seed)

        # Note that no shuffling of order of genome occurs. This makes the problem "harder" because for the most
        # part we're training / testing on non-adjacent regions (except at the boundaries)
        cv_results = cross_validate(model, X, y, cv=10, return_estimator=True, scoring="r2", n_jobs=-1)
        mean_score = np.mean(cv_results["test_score"])
        if mean_score > self._best_score:
            self._best_score = mean_score
            self.best_model_perfoldtrained = cv_results["estimator"]
            print("Saving new best optimizer...")
            pickle.dump({"best_score": self._best_score,
                         "best_model_perfoldtrained": self.best_model_perfoldtrained},
                        open(self._savepath, "wb"))
            print("Done saving")
        return mean_score


    def optimize_optuna_study(self, study_name, ML_model, X_train, y_train, seed, n_optuna_trials,
                              sqlite):
        storage_name = get_storage_name(sqlite)

        study = optuna.create_study(direction="maximize",
                                    storage=storage_name,
                                    study_name=study_name,
                                    load_if_exists=True,
                                    sampler=optuna.samplers.TPESampler(seed=seed))
        n_existing_complete_trials = len([trial for trial in study.trials if trial.state ==
                                            optuna.trial.TrialState.COMPLETE])
        print(f"Number of existing COMPLETE optuna trials: {n_existing_complete_trials}")
        n_optuna_trials_remaining = n_optuna_trials - n_existing_complete_trials
        n_optuna_trials_remaining = max(0, n_optuna_trials_remaining)
        if n_optuna_trials_remaining > 0:
            print(f"Running an extra {n_optuna_trials_remaining} trials")
            study.optimize(lambda trial: self._optuna_objective(trial, ML_model=ML_model, X=X_train,
                                                                y=y_train, seed=seed),
                           n_trials=n_optuna_trials_remaining)
        else:
            print(f"Done running {n_optuna_trials} trials!")
        return study


def get_and_save_test_set_perf(X_test, y_test, model, filepath):
    test_preds = model.predict(X_test)
    test_set_performance = r2_score(y_test, test_preds)
    with open(filepath, "w") as f:
        f.write(str(test_set_performance))
    return test_set_performance


def train_val_test(scATAC_df, mutations, backwards_elim_dir, test_set_perf_filepath,
                   ML_model, seed, scATAC_dir, cancer_type_or_donor_id,
                   n_optuna_trials_prebackward_selection, n_optuna_trials_backward_selection,
                   feature_importance_method, sqlite, debug_bfs, fold_for_test_set):
    X_train, X_test, y_train, y_test = get_train_test_split(scATAC_df, mutations, 10, fold_for_test_set)
    if seed == 1:
        X_train.to_csv(f"{os.path.dirname(backwards_elim_dir)}/X_train.csv")
        X_test.to_csv(f"{os.path.dirname(backwards_elim_dir)}/X_test.csv")
        y_train.to_csv(f"{os.path.dirname(backwards_elim_dir)}/y_train.csv")
        y_test.to_csv(f"{os.path.dirname(backwards_elim_dir)}/y_test.csv")

    # Define as None in case scATAC_df.shape[1] <= 20
    n = None
    best_model_fulldatatrained = None
    best_model_perfoldtrained = None

    if scATAC_df.shape[1] > 20:
        print("Getting a starter model!")
        study_name = f"{cancer_type_or_donor_id}_{scATAC_dir}"
        model_optimizer = ModelOptimizer(backwards_elim_dir + "/" + "model_optimizer_starter.pkl")
        study = model_optimizer.optimize_optuna_study(study_name=study_name, ML_model=ML_model, X_train=X_train,
                                                      y_train=y_train, seed=seed,
                                                      n_optuna_trials=n_optuna_trials_prebackward_selection,
                                                      sqlite=sqlite)
        print("Done getting starter model!")
        best_params = study.best_params
        if os.path.exists(f"{backwards_elim_dir}/best_model_fulldatatrained.pkl"):
            print("Loading existing full data trained model...")
            best_model_fulldatatrained = pickle.load(open(f"{backwards_elim_dir}/best_model_fulldatatrained.pkl", "rb"))
        else:
            if ML_model == "XGB":
                best_model_fulldatatrained = XGBRegressor(**best_params)
            print("Training full train data model...")
            best_model_fulldatatrained.fit(X=X_train, y=y_train)
            pickle.dump(best_model_fulldatatrained, open(f"{backwards_elim_dir}/best_model_fulldatatrained.pkl", "wb"))

        best_model_perfoldtrained = model_optimizer.best_model_perfoldtrained
        # Test Set Performance
        print("Getting test set performance using full feature space...")
        tsp = get_and_save_test_set_perf(X_test, y_test, best_model_fulldatatrained, test_set_perf_filepath)
        print(f"Test set performance with all features: {tsp}")
        # For backward feature selection
        n = 20
        filepath = f"top_features_iteration_{n - 1}"
    else:
        print("Starter model not needed! Number of features is less than or equal to 20 already!")
        print(X_train.shape[1])
        filepath = f"top_features_iteration_{X_train.shape[1] - 1}"

    if feature_importance_method != "default_importance":
        filepath = filepath + f"_by_{feature_importance_method}"

    filepath = filepath + ".txt"

    if not os.path.exists(f"{backwards_elim_dir}/{filepath}") or debug_bfs:
        print("Running backward feature selection...")
        backward_eliminate_features(X_train, y_train, backwards_elim_dir, ML_model, scATAC_dir,
                                    cancer_type_or_donor_id, seed, n_optuna_trials_backward_selection,
                                    feature_importance_method,
                                    sqlite,
                                    best_model_fulldatatrained=best_model_fulldatatrained,
                                    best_model_perfoldtrained=best_model_perfoldtrained,
                                    starting_n=n)
    else:
        print("Backward feature selection is already done!")


def save_model_with_n_features_test_performance(scATAC_df, mutations_df, scATAC_dir, n, ML_model, cancer_type,
                                                feature_importance_method, fold_for_test_set):
    model = load_n_features_backwards_elim_models(n, scATAC_df.shape[1], cancer_type, ML_model, scATAC_dir,
                                                  feature_importance_method)
    # print(scATAC_df.shape[1])
    if scATAC_df.shape[1] > 20:
       model_iteration = 20 - n + 1
    else:
       model_iteration = scATAC_df.shape[1] - n + 1

    scATAC_df = scATAC_df.loc[:, model.feature_names_in_]
    scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
    _, X_test, _, y_test = get_train_test_split(scATAC_df, mutations_df, 10, fold_for_test_set)
    test_set_perf_filepath = f"models/{ML_model}/" \
                             f"{cancer_type}/{scATAC_dir}/backwards_elimination_results/" \
                             f"model_iteration_{model_iteration}"
    if feature_importance_method != "default_importance":
        test_set_perf_filepath = test_set_perf_filepath + f"_feature_importance_{feature_importance_method}"

    test_set_perf_filepath = test_set_perf_filepath + "_test_performance.txt"
    tsp = get_and_save_test_set_perf(X_test, y_test, model, test_set_perf_filepath)
    print(f"Test set performance with {n} features: {tsp}")


def apply_func_to_kfolds(X, y, func, **kwargs):
    kf = KFold(n_splits=10)
    func_out = []
    for i, (train_index, val_index) in enumerate(kf.split(X)):
        X_train, y_train = X.iloc[train_index], y.iloc[train_index]
        X_valid, y_valid = X.iloc[val_index], y.iloc[val_index]
        kwargs["i"] = i
        kwargs["X_train"] = X_train
        kwargs["y_train"] = y_train
        func_out.append(func(X_valid, y_valid, **kwargs))
    return func_out

#### Call other scripts ####
# def call_plot_top_features(seed, cancer_types_arg, ML_model, datasets_arg, scATAC_cell_number_filter,
#                            annotation_dir, top_features_to_plot, feature_importance_method, fold_for_test_set):
#     print(f"Plotting top features for seed {seed}...")
#     command = ["Rscript", "plot_top_features.R",
#                          f"--cancer_types={cancer_types_arg}",
#                          f"--ML_model={ML_model}",
#                          f"--datasets={datasets_arg}",
#                          f"--seed={seed}",
#                          f"--cell_number_filter={scATAC_cell_number_filter}",
#                          f"--annotation={annotation_dir}",
#                          f"--top_features_to_plot={','.join(list(map(str, top_features_to_plot)))}",
#                          f"--feature_importance_method={feature_importance_method}",
#                          f"--fold_for_test_set={fold_for_test_set}"]
#     subprocess.call(command)
#     print(f"Done plotting top features for seed {seed}!")

def call_plot_top_features(seed_range, cancer_types_arg, ML_model, datasets_arg, scATAC_cell_number_filter,
                           annotation_dir, top_features_to_plot, feature_importance_method, fold_for_test_set):
    print(f"Plotting top features for seed range {seed_range}...")
    command = ["Rscript", "plot_top_features.R",
                         f"--cancer_types={cancer_types_arg}",
                         f"--ML_model={ML_model}",
                         f"--datasets={datasets_arg}",
                         f"--seed_range={'-'.join(list(map(str, top_features_to_plot)))}",
                         f"--cell_number_filter={scATAC_cell_number_filter}",
                         f"--annotation={annotation_dir}",
                         f"--top_features_to_plot={','.join(list(map(str, top_features_to_plot)))}",
                         f"--feature_importance_method={feature_importance_method}",
                         f"--folds_for_test_set={'-'.join(list(map(str, [fold_for_test_set, fold_for_test_set])))}"]
    subprocess.call(command)
    print(f"Done plotting top features for seed range {seed_range}!")

#### Other ####
def construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir, seed,
                         fold_for_test_set):
    scATAC_dir = f"scATAC_source_{scATAC_sources}_cell_number_filter_{scATAC_cell_number_filter}"
    if tss_filter:
        scATAC_dir = scATAC_dir + "_tss_fragment_filter_" + tss_filter
    scATAC_dir = scATAC_dir + f"_annotation_{annotation_dir}_seed_{seed}_fold_for_test_set_{fold_for_test_set + 1}"
    return scATAC_dir


def construct_scATAC_sources(datasets):
    scATAC_sources = ""

    for idx, dataset in enumerate(datasets):
        if scATAC_sources == "":
            scATAC_sources = dataset
        else:
            scATAC_sources = "_".join((scATAC_sources, dataset))
    return scATAC_sources


def get_storage_name(sqlite=False):
    if sqlite:
        storage_name = "sqlite:///temp.db"
    else:
        hostname_file = open(os.path.dirname(os.path.abspath(__file__)) + "/" + "postgresql_hostname.txt", "r")
        hostname = hostname_file.readline().strip()
        storage_name = f"postgresql://bgiotti:bgiotti@{hostname}:5432/optuna_db"
    return storage_name


def load_n_features_backwards_elim_models(n, total_num_features, cancer_type, ML_model, scATAC_dir,
                                         feature_importance_method, full_data_trained=True):
    print(n)
    if total_num_features > 20:
        model_iteration = 20 - n + 1
    else:
        model_iteration = total_num_features - n + 1
    if full_data_trained:
        filename = f"model_iteration_{model_iteration}"
    else:
        filename = f"per_fold_model_iteration_{model_iteration}"

    if feature_importance_method != "default_importance":
        filename = filename + f"_feature_importance_{feature_importance_method}"
    filename = filename + ".pkl"

    backwards_elim_model_file = f"models/{ML_model}/" \
                                f"{cancer_type}/{scATAC_dir}/backwards_elimination_results/{filename}"
    model = pickle.load(open(backwards_elim_model_file, "rb"))
    return model

def compute_error(X_val, y_val, **kwargs):
    estimator = kwargs["estimator"]
    if kwargs["train_new_model"]:
        estimator = XGBRegressor(**estimator.get_params())
        estimator.fit(X=kwargs["X_train"], y=kwargs["y_train"])
    if kwargs["i"] is not None:
        print(f"Calculating error for fold {kwargs['i'] + 1}...")
    preds = estimator.predict(X_val)
    abs_err = np.absolute(y_val - preds)
    percent_error = abs_err / (y_val + 1) * 100
    return {"abs_err": abs_err, "percent_err": percent_error, "preds": preds}