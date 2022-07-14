#### Imports ####
import pyreadr
import pandas as pd
import numpy as np
import pickle
import os
import time
from pipelinehelper import PipelineHelper
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.metrics import explained_variance_score
import argparse
import glob

parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', required=True)
parser.add_argument('--all_cells', action="store_true",
                    help='run model on all cells', default=False)
parser.add_argument('--tissue_spec_cells', action="store_true",
                    help='run model on tissue specific cells', default=False)
parser.add_argument('--clustered_mutations', action="store_true",
                    help='run model on hierarchically clustered mutations', default=False)
parser.add_argument('--bing_ren', action="store_true",
                    help='use Bing Ren ATACseq', default=False)
parser.add_argument('--shendure', action="store_true",
                    help='use Shendure ATACseq', default=False)

config = parser.parse_args()
cancer_types = config.cancer_types
run_all_cells = config.all_cells
run_tissue_spec_cells = config.tissue_spec_cells
run_clustered_mutations = config.clustered_mutations
bing_ren = config.bing_ren
shendure = config.shendure

#### Helpers ####
# Load Data helpers
def load_scATAC(scATAC_path):
    scATAC_df = pyreadr.read_r(scATAC_path)
    scATAC_df = scATAC_df[None]
    scATAC_df = scATAC_df.T
    return scATAC_df

def load_agg_mutations():
    return pd.read_csv("processed_data/mut_count_data.csv",
                       index_col=0)

# Filter Data helpers
def filter_agg_data(scATAC_df, mutations_df):
    idx_select = ~pd.isnull(mutations_df).any(axis=1)
    scATAC_df = scATAC_df[idx_select]
    mutations_df = mutations_df[idx_select]
    return scATAC_df, mutations_df

def filter_clustered_data(scATAC_df, mutations_df):
    return pd.DataFrame(scATAC_df, index=mutations_df.index)
# Split Train/Test helpers
def get_train_test_split(X, y, test_size):
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=test_size, random_state=42)
    return X_train, X_test, y_train, y_test

# Cross Validation helpers
def grid_search(X_train, y_train, pipe, params, num_k_folds):
    start_time = time.time()
    grid_search_object = GridSearchCV(pipe, params, scoring="explained_variance",
                                      cv=KFold(num_k_folds), n_jobs=-1, verbose=100)
    grid_search_object.fit(X_train, y_train)
    print(f"--- {time.time() - start_time} seconds ---")
    return grid_search_object

# Feature Selection helpers
def get_top_n_features(clf, n, features):
    feature_importances = clf.feature_importances_
    feat_importance_idx = np.argsort(feature_importances)[::-1]
    top_n_feats = features[feat_importance_idx][:n]
    return top_n_feats

def print_and_save_top_features(top_features, filepath):
    n = len(top_features)
    print(f"Top {n} features")
    f = open(filepath, "w")
    for idx, feature in enumerate(top_features):
        print(f"{idx+1}. {feature}")
        f.write(f"{idx+1}. {feature}\n")

def backward_eliminate_features(X_train, y_train, starting_clf, starting_n,
                                params, num_k_folds, backwards_elim_dir):
    os.makedirs(backwards_elim_dir, exist_ok=True)
    top_n_feats = get_top_n_features(starting_clf, starting_n, X_train.columns.values)
    X_train = X_train.loc[:, top_n_feats]
    print_and_save_top_features(top_n_feats,
                                filepath=f"{backwards_elim_dir}/starter_model_top_features.txt")
    for idx in range(1, len(top_n_feats)):
        filepath = f"{backwards_elim_dir}/top_features_iteration_{idx}.txt"
        if (not os.path.exists(filepath)):
            pipe = Pipeline([
                ('regressor', PipelineHelper([
                    ('rf', RandomForestRegressor(random_state=idx)),
                ])),
            ])
            grid_search_results = grid_search(X_train, y_train, pipe, params, num_k_folds)
            pickle.dump(grid_search_results, open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'wb'))
        grid_search_results = pickle.load(open(f"{backwards_elim_dir}/model_iteration_{idx}.pkl", 'rb'))
        best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
        top_n_feats = get_top_n_features(best_model, len(X_train.columns) - 1, X_train.columns)
        X_train = X_train.loc[:, top_n_feats]
        if (not os.path.exists(filepath)):
            print_and_save_top_features(top_n_feats, filepath=filepath)


# Test set performance helpers
def print_and_save_test_set_perf(X_test, y_test, model, filename):
    test_preds = model.predict(X_test)
    test_set_performance = explained_variance_score(y_test, test_preds)
    with open(filename, "w") as f:
        f.write(str(test_set_performance))
    print(f"Test set performance: {test_set_performance}")

def train_val_test(scATAC_df, mutations, cv_filename, backwards_elim_dir,
                   test_set_perf_filename):
    X_train, X_test, y_train, y_test = get_train_test_split(scATAC_df, mutations,
                                                            0.10)
    pipe = Pipeline([
        ('regressor', PipelineHelper([
            ('rf', RandomForestRegressor(random_state=0)),
        ])),
    ])

    params = {
        'regressor__selected_model': pipe.named_steps['regressor'].generate({
            'rf__n_estimators':[10, 100, 1000],
        })
    }

    if (not os.path.exists(cv_filename)):
        grid_search_results = grid_search(X_train, y_train, pipe, params, 10)
        pickle.dump(grid_search_results, open(cv_filename, 'wb'))

    grid_search_results = pickle.load(open(cv_filename, 'rb'))
    print(f"Best Score: {grid_search_results.best_score_}\n")
    best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
    n = 20
    # top_n_feats = get_top_n_features(best_model, n, X_train.columns.values)
    # top_n_feats_tissue_spec = get_top_n_features(best_model_tissue_spec, n, X_train_tissue_spec.columns.values)
    #print_top_features(top_n_feats)

    backward_eliminate_features(X_train, y_train, best_model, n, params, 10, backwards_elim_dir)

    #### Test Set Performance ####
    print_and_save_test_set_perf(X_test, y_test, best_model, test_set_perf_filename)

def run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec, scATAC_source="bing_ren"):
    mutations_df = load_agg_mutations()

    #### Filter data ####
    scATAC_df, mutations_df = filter_agg_data(scATAC_df, mutations_df)

    mutations_per_cancer_type = {}
    for cancer_type in cancer_types:
        os.makedirs(f"models/{cancer_type}/scATAC_source_{scATAC_source}", exist_ok=True)
        cancer_type_specific_mut_idx = [i for i, s in enumerate(mutations_df.columns.values) if cancer_type == s][0]
        mutations_per_cancer_type[cancer_type] = mutations_df.iloc[:, cancer_type_specific_mut_idx]

        if (run_all_cells):
            backwards_elim_dir=f"models/{cancer_type}/scATAC_source_{scATAC_source}/backwards_elimination_results"
            grid_search_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/grid_search_results.pkl"
            test_set_perf_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/test_set_performance.txt"

            # All Cells
            train_val_test(scATAC_df, mutations_per_cancer_type[cancer_type],
                           grid_search_filename,
                           backwards_elim_dir,
                           test_set_perf_filename)

        # Tissue Specific
        if (run_tissue_spec):
            backwards_elim_dir=f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
                               f"backwards_elimination_results_tissue_spec"
            grid_search_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
                                   f"grid_search_results_tissue_specific.pkl"
            test_set_perf_filename = f"models/{cancer_type}/scATAC_source_{scATAC_source}/" \
                                     f"test_set_performance_tissue_spec.txt"

            tissue = cancer_type.split("-")[0]
            tissue_specific_cell_types = [cell_type for cell_type in scATAC_df.columns.values if tissue in cell_type]
            per_tissue_df = scATAC_df.loc[:, tissue_specific_cell_types]
            train_val_test(per_tissue_df, mutations_per_cancer_type[cancer_type],
                           grid_search_filename,
                           backwards_elim_dir,
                           test_set_perf_filename)


def run_clustered_data_analysis(scATAC_df):
    for cancer_type in cancer_types:
        os.makedirs(f"models/{cancer_type}", exist_ok=True)
        cancer_hierarchical_dir = f"processed_data/hierarchically_clustered_mutations/{cancer_type}"
        for cluster_method_dir in os.listdir(cancer_hierarchical_dir):
            for threshold_dir in os.listdir(f"{cancer_hierarchical_dir}/{cluster_method_dir}"):
                run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
                                       threshold_dir)

def run_per_cluster_models(scATAC_df, cancer_type, cancer_hierarchical_dir, cluster_method_dir,
                           threshold_dir):
    for cluster_file in glob.glob(f"{cancer_hierarchical_dir}/{cluster_method_dir}/{threshold_dir}/aggregated*"):
        mutations_df = pd.read_csv(cluster_file, index_col=0)
        cluster_dir = cluster_file.split("/")[-1].split(".")[0]
        cluster_dir = f"models/{cancer_type}/{cluster_method_dir}/{threshold_dir}/{cluster_dir}/"
        os.makedirs(cluster_dir, exist_ok=True)
        backwards_elim_dir=f"{cluster_dir}/backwards_elimination_results"
        grid_search_filename = f"{cluster_dir}/grid_search_results.pkl"
        test_set_perf_filename = f"{cluster_dir}/test_set_performance.txt"

        scATAC_df = filter_clustered_data(scATAC_df, mutations_df)
        mutations = mutations_df.values.reshape(-1)
        train_val_test(scATAC_df, mutations,
                       grid_search_filename,
                       backwards_elim_dir,
                       test_set_perf_filename)

#### Load scATAC ####
scATAC_df = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
                        "/count_filter_100_combined_count_overlaps.rds")
if ((run_all_cells or run_tissue_spec_cells) and bing_ren):
    run_unclustered_data_analysis(scATAC_df, run_all_cells, run_tissue_spec_cells)
if (run_clustered_mutations and bing_ren):
    run_clustered_data_analysis(scATAC_df)

scATAC_df_shendure = load_scATAC("processed_data/count_overlap_data/combined_count_overlaps" \
                                 "/shendure_count_filter_100_combined_count_overlaps.rds")
if ((run_all_cells or run_tissue_spec_cells) and shendure):
    run_unclustered_data_analysis(scATAC_df_shendure, run_all_cells, run_tissue_spec_cells, "shendure")