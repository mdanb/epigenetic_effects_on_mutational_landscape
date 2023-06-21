import pyreadr
from natsort import natsorted
from itertools import chain
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import pickle
import os
from sklearn.metrics import r2_score
from xgboost import XGBRegressor
# from sklearn.model_selection import cross_val_score
import optuna
import xgboost as xgb


### Load Data helpers ###
def load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                   histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                   CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor):
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
        mutations_df = load_de_novo_seurat_clustered_cancers(cancer_types)
    elif CPTAC:
        mutations_df = load_CPTAC()
    elif combined_CPTAC_ICGC:
        mutations_df = load_combined_CPTAC_ICGC()
    elif RNA_subtyped:
        mutations_df = load_RNA_subtyped_mutations()
    elif per_donor:
        mutations_df = load_per_donor_mutations(cancer_types[0])
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
    return(df.loc[natsorted(df.index)])

def load_CPTAC():
    df = pd.read_csv(f"../../data/processed_data/CPTAC_mutations.csv",
                     index_col=0)
    return(df.loc[natsorted(df.index)])

def load_combined_CPTAC_ICGC():
    df = pd.read_csv(f"../../data/processed_data/combined_CPTAC_ICGC_mutations.csv",
                     index_col=0)
    return(df.loc[natsorted(df.index)])

def load_agg_mutations():
    df = pd.read_csv("../../data/processed_data/mut_count_data.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_woo_pcawg_mutations():
    df = pd.read_csv("../../data/processed_data/pcawg_agg_woo.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_RNA_subtyped_mutations():
    df = pd.read_csv("../../data/processed_data/RNA_subtyped_cancers.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_meso():
    df = pd.read_csv("../../data/processed_data/mesothelioma.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_sclc_mutations():
    df = pd.read_csv("../../data/processed_data/sclc_count_overlaps.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_subtyped_lung_mutations():
    df = pd.read_csv("../../data/processed_data/lung_subtyped_mutation.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_histologically_subtyped_mutations():
    df = pd.read_csv("../../data/processed_data/histologically_subtyped_mutations.csv",
                       index_col=0)
    return(df.loc[natsorted(df.index)])

def load_de_novo_seurat_clustered_cancers(cancer_types):
    print("Loading De Novo Seurat clustered cancers")
    cancer_type = cancer_types[0].split("x")[0]
    seurat_cluster_settings = cancer_types[0].split("x")[1]
    df = pd.read_csv(f"../../data/processed_data/de_novo_seurat_clustered_mutations/{cancer_type}/"
                     f"{seurat_cluster_settings}.csv", index_col=0)
    return(df.loc[natsorted(df.index)])

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
    return(scATAC_df)

#### Dataframe curation helpers ####
def add_na_ranges(mutations_df):
    full_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
    mutations_df = pd.merge(full_ranges, mutations_df, left_on = "x", right_index=True,
                            how="outer").set_index("x")
    return(mutations_df.loc[natsorted(mutations_df.index)])

def add_dataset_origin_to_cell_types(df, dataset):
    if (dataset == "Bingren"):
        df.columns = [c + " BR" for c in df.columns]
    elif (dataset == "Shendure"):
        df.columns = [c + " SH" for c in df.columns]
    elif (dataset == "Tsankov"):
        df.columns = [c + " TS" for c in df.columns]
    elif (dataset == "Greenleaf_brain"):
        df.columns = [c + " GL_Br" for c in df.columns]
    elif (dataset == "Greenleaf_pbmc_bm"):
        df.columns = [c + " GL_BlBm" for c in df.columns]
    elif (dataset == "Yang_kidney"):
        df.columns = [c + " Y_K" for c in df.columns]
    return(df)

def construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir):
    datasets_combined_count_overlaps = []
    for dataset in datasets:
        if (tss_filter):
            tss_filtered_root = "../../data/processed_data/count_overlap_data/tsse_filtered"
            chr_ranges = pd.read_csv("../../data/processed_data/chr_ranges.csv")
            scATAC_df = load_scATAC(f"{tss_filtered_root}/{dataset}/combined/{annotation_dir}/" \
                                    f"combined_{tss_filter}_fragments.rds").T
            scATAC_df.index = chr_ranges["x"].values
            datasets_combined_count_overlaps.append(scATAC_df)
        else:
            scATAC_df = load_scATAC("../../data/processed_data/count_overlap_data/combined_count_overlaps" \
            f"/{annotation_dir}/{dataset}_combined_count_overlaps.rds")
            metadata = load_scATAC_metadata("../../data/processed_data/count_overlap_data/combined_count_overlaps" \
            f"/{annotation_dir}/{dataset}_combined_count_overlaps_metadata.rds")
            scATAC_df = filter_scATAC_df_by_num_cell_per_cell_type(scATAC_df, scATAC_cell_number_filter, metadata)
            datasets_combined_count_overlaps.append(scATAC_df)

    for idx, dataset in enumerate(datasets):
        datasets_combined_count_overlaps[idx] = [add_dataset_origin_to_cell_types(datasets_combined_count_overlaps[idx],
                                                 dataset)]
    scATAC_df = pd.concat(chain(*datasets_combined_count_overlaps), axis=1)
    return(scATAC_df)


#### Split Train/Test helpers ####
def get_train_test_split(X, y, test_size, seed):
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=test_size, random_state=seed)
    return X_train, X_test, y_train, y_test

#### Feature Selection helpers ####
def get_top_n_features(clf, n, features):
    feature_importances = clf.feature_importances_
    feat_importance_idx = np.argsort(feature_importances)[::-1]
    top_n_feats = features[feat_importance_idx][:n]
    return top_n_feats

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
            f.write(f"-{feature}\n")

def backward_eliminate_features(X_train, y_train, backwards_elim_dir,
                                ML_model, scATAC_dir, cancer_type_or_donor_id, seed,
                                n_optuna_trials_backward_selection, starting_clf=None, starting_n=None):
    os.makedirs(backwards_elim_dir, exist_ok=True)
    if X_train.shape[1] > 20:
        top_n_feats = get_top_n_features(starting_clf, starting_n, X_train.columns.values)
        all_feature_rankings = get_top_n_features(starting_clf, len(X_train.columns.values), X_train.columns.values)
        print_and_save_features(all_feature_rankings, filepath=f"{backwards_elim_dir}/all_features_rankings.txt",
                                top=True)
        X_train = X_train.loc[:, top_n_feats]
        print_and_save_features(top_n_feats,
                                filepath=f"{backwards_elim_dir}/starter_model_top_features.txt",
                                top=True)
        num_iterations = len(top_n_feats)
    else:
        print_and_save_features(X_train.columns,
                                filepath=f"{backwards_elim_dir}/starter_features.txt",
                                top=False)
        num_iterations = X_train.shape[1]

    for idx in range(1, num_iterations):
        filepath = f"{backwards_elim_dir}/top_features_iteration_{idx}.txt"
        # if not os.path.exists(filepath):
        study = optimize_optuna_study(study_name=f"{cancer_type_or_donor_id}_iter_{idx}_{scATAC_dir}",
                                      ML_model=ML_model, X_train=X_train,
                                      y_train=y_train,
                                      seed=seed,
                                      n_optuna_trials=n_optuna_trials_backward_selection)
        best_params = study.best_params
        if ML_model == "XGB":
            best_model = XGBRegressor(**best_params)

        best_model.fit(X=X_train, y=y_train)
        pickle.dump(best_model, open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'wb'))
        # grid_search_results = pickle.load(open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'rb'))
        top_n_feats = get_top_n_features(best_model, len(X_train.columns) - 1, X_train.columns)
        X_train = X_train.loc[:, top_n_feats]
        if not os.path.exists(filepath):
            print_and_save_features(top_n_feats, filepath=filepath, top=True)

#### Model train/val/test helpers ####
def optimize_optuna_study(study_name, ML_model, X_train, y_train, seed, n_optuna_trials):
    # storage_name = "mysql+pymysql://mdanb:mdanb@localhost:3306/optuna_db"
    #storage_name = "mysql+pymysql://mdanb:mdanb@localhost:3306/optuna_db"
    storage_name = "sqlite:///example.db"
    study = optuna.create_study(direction="maximize",
                                storage=storage_name,
                                study_name=study_name,
                                load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=seed),
                                pruner=optuna.pruners.MedianPruner(n_warmup_steps=5))
    n_existing_trials = len(study.trials)
    print(f"Number of existing optuna trials: {n_existing_trials}")
    n_optuna_trials = n_optuna_trials - n_existing_trials
    n_optuna_trials_remaining = max(0, n_optuna_trials)
    if n_optuna_trials > 0:
        print(f"Running an extra {n_optuna_trials_remaining} trials")
    else:
        print(f"Done running {n_optuna_trials} trials!")
    study.optimize(lambda trial: optuna_objective(trial, ML_model=ML_model, X=X_train, y=y_train,
                                                  seed=seed), n_trials=n_optuna_trials_remaining)
    return study


    # scores = []
    # kf = KFold(n_splits=10, shuffle=True, random_state=seed)
    #
    # if ML_model == "XGB":
    #     param = {
    #         'max_depth': trial.suggest_int('max_depth', 3, 10),
    #         'learning_rate': trial.suggest_float('learning_rate', 1e-8, 1.0, log=True),
    #         'subsample': trial.suggest_float('subsample', 0.1, 1.0),
    #         'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
    #         'min_child_weight': trial.suggest_int('min_child_weight', 1, 6),
    #         'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
    #         'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
    #     }
    #     num_boost_round = trial.suggest_int('num_boost_round', 100, 500)
    #
    #     for train_index, val_index in kf.split(X):
    #         X_train, X_val = X.iloc[train_index], X.iloc[val_index]
    #         y_train, y_val = y[train_index], y[val_index]
    #
    #         dtrain = xgb.DMatrix(X_train, label=y_train)
    #         dval = xgb.DMatrix(X_val, label=y_val)
    #
    #         watchlist = [(dtrain, 'train'), (dval, 'eval')]
    #
    #         # model = xgb.train(param, dtrain, num_boost_round=num_boost_round, evals=watchlist,
    #         #                   callbacks=[optuna.integration.XGBoostPruningCallback(trial, "eval-rmse")],
    #         #                   early_stopping_rounds=10, verbose_eval=False)
    #         model = xgb.train(param, dtrain, num_boost_round=num_boost_round, evals=watchlist,
    #                           early_stopping_rounds=10, verbose_eval=False)
    #
    #         preds = model.predict(dval)
    #         score = r2_score(y_val, preds)
    #         scores.append(score)
    # return np.mean(scores)


def optuna_objective(trial, ML_model, X, y, seed):
    if ML_model == "XGB":
        param = {
            'max_depth': trial.suggest_int('max_depth', 3, 10),
            'learning_rate': trial.suggest_float('learning_rate', 1e-8, 1.0, log=True),
            'subsample': trial.suggest_float('subsample', 0.1, 1.0),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
            'min_child_weight': trial.suggest_int('min_child_weight', 1, 6),
            'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
            'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
            'objective': 'reg:squarederror',  # Replace with your objective
            'seed': seed,
        }
        num_boost_round = trial.suggest_int('num_boost_round', 100, 500)

        dtrain = xgb.DMatrix(X, label=y)

        cv_results = xgb.cv(param, dtrain, num_boost_round=num_boost_round,
                        nfold=10, stratified=False,
                        early_stopping_rounds=10,
                        callbacks=[optuna.integration.XGBoostPruningCallback(trial, "test-rmse")],
                        seed=seed)

    return cv_results['test-rmse-mean'].values[-1]

def print_and_save_test_set_perf(X_test, y_test, model, filepath):
    test_preds = model.predict(X_test)
    test_set_performance = r2_score(y_test, test_preds)
    with open(filepath, "w") as f:
        f.write(str(test_set_performance))
    print(f"Test set performance: {test_set_performance}")

def train_val_test(scATAC_df, mutations, backwards_elim_dir, test_set_perf_filepath,
                   ML_model, seed, scATAC_dir, cancer_type_or_donor_id,
                   n_optuna_trials_prebackward_selection, n_optuna_trials_backward_selection):
    X_train, X_test, y_train, y_test = get_train_test_split(scATAC_df, mutations, 0.10, seed)

    starting_clf = None
    n = None

    if scATAC_df.shape[1] > 20:
        print("Getting a starter model!")
        study_name = f"{cancer_type_or_donor_id}_{scATAC_dir}"
        study = optimize_optuna_study(study_name=study_name, ML_model=ML_model, X_train=X_train, y_train=y_train, seed=seed,
                                      n_optuna_trials=n_optuna_trials_prebackward_selection)
        best_params = study.best_params
        if ML_model == "XGB":
            starting_clf = XGBRegressor(**best_params)
        n = 20
        starting_clf.fit(X=X_train, y=y_train)
        # Test Set Performance
        print_and_save_test_set_perf(X_test, y_test, starting_clf, test_set_perf_filepath)
        print("Done getting starter model!")
    else:
        print("Starter model not needed! Number of features is less than or equal to 20 already!")
    if not os.path.exists(f"{backwards_elim_dir}/top_features_iteration_{scATAC_df.shape[1] - 1}.txt"):
        backward_eliminate_features(X_train, y_train, backwards_elim_dir, ML_model, scATAC_dir,
                                    cancer_type_or_donor_id, seed, n_optuna_trials_backward_selection,
                                    starting_clf=starting_clf, starting_n=n)


def save_n_features_model_test_performance(n, datasets, ML_model, scATAC_cell_number_filter, tss_filter, annotation_dir,
                                           meso, SCLC, lung_subtyped, woo_pcawg,
                                           histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                                           CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor, seed):
    for cancer_type in cancer_types:
        scATAC_sources = construct_scATAC_sources(datasets)
        scATAC_dir = construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter,
                                         annotation_dir, seed)
        filename = f"model_iteration_{n}.pkl"
        backwards_elim_model_file = f"models/{ML_model}/" \
                                    f"{cancer_type}/{scATAC_dir}/backwards_elimination_results/{filename}"
        model = pickle.load(open(backwards_elim_model_file, "rb"))
        scATAC_df = construct_scATAC_df(tss_filter, datasets, scATAC_cell_number_filter, annotation_dir)
        scATAC_df = scATAC_df.loc[:, model.feature_names_in_]
        scATAC_df = scATAC_df.loc[natsorted(scATAC_df.index)]
        mutations_df = load_mutations(meso, SCLC, lung_subtyped, woo_pcawg,
                                      histologically_subtyped_mutations, de_novo_seurat_clustering, cancer_types,
                                      CPTAC, combined_CPTAC_ICGC, RNA_subtyped, per_donor)

        if not pd.isna(mutations_df).any().any():
            # for compatibility
            mutations_df = add_na_ranges(mutations_df)
        scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)
        cancer_specific_mutations = filter_mutations_by_cancer(mutations_df, cancer_type)

        _, X_test, _, y_test = get_train_test_split(scATAC_df, cancer_specific_mutations, 0.10, seed)
        test_set_perf_filepath = f"models/{ML_model}/" \
                                 f"{cancer_type}/{scATAC_dir}/backwards_elimination_results/" \
                                 f"model_iteration_{n}_test_performance.txt"
        print_and_save_test_set_perf(X_test, y_test, model, test_set_perf_filepath)

#### Other ####
def construct_scATAC_dir(scATAC_sources, scATAC_cell_number_filter, tss_filter, annotation_dir, seed):
    scATAC_dir = f"scATAC_source_{scATAC_sources}_cell_number_filter_{scATAC_cell_number_filter}"
    if (tss_filter):
        scATAC_dir = scATAC_dir + "_tss_fragment_filter_" + tss_filter
    scATAC_dir = scATAC_dir + f"_annotation_{annotation_dir}_seed_{seed}"
    return(scATAC_dir)

def construct_scATAC_sources(datasets):
    scATAC_sources = ""

    for idx, dataset in enumerate(datasets):
        if (scATAC_sources == ""):
            scATAC_sources = dataset
        else:
            scATAC_sources = "_".join((scATAC_sources, dataset))
    return(scATAC_sources)
