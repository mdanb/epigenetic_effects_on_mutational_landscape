# Imports
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
parser = argparse.ArgumentParser()

parser.add_argument('--cancer_types', nargs="+", type=str,
                    help='which cancer types to analyze', required=True)

config = parser.parse_args()
cancer_types = config.cancer_types

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
                                params, num_k_folds, cancer_type):
    dir=f"models/{cancer_type}/backwards_elimination_results"
    os.makedirs(dir, exist_ok=True)
    top_n_feats = get_top_n_features(starting_clf, starting_n, X_train.columns.values)
    X_train = X_train.loc[:, top_n_feats]
    print_and_save_top_features(top_n_feats,
                                filepath=f"{dir}/starter_model_top_features.txt")
    idx = 1
    while (len(top_n_feats) > 1):
        pipe = Pipeline([
            ('regressor', PipelineHelper([
                ('rf', RandomForestRegressor(random_state=idx)),
            ])),
        ])
        grid_search_results = grid_search(X_train, y_train, pipe, params, num_k_folds)
        best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
        top_n_feats = get_top_n_features(best_model, len(X_train.columns) - 1, X_train.columns)
        X_train = X_train.loc[:, top_n_feats]
        print_and_save_top_features(top_n_feats, filepath=f"{dir}/top_features_iteration_{idx}.txt")
        idx += 1

# Load data
scATAC_df = pyreadr.read_r("count_overlap_data/processed_count_overlaps" \
                           "/count_filter_100_combined_count_overlaps.rds")
scATAC_df = scATAC_df[None]
scATAC_df = scATAC_df.T
mutations_df = pd.read_csv("processed_data/mut_count_data.csv", index_col=0)

# Select data
idx_select = ~pd.isnull(mutations_df).any(axis=1)
scATAC_df = scATAC_df[idx_select]
mutations_df = mutations_df[idx_select]

mutations_per_cancer_type = {}
for cancer_type in cancer_types:
    # IF FILES DON'T EXIST
    cancer_type_specific_mut_idx = [i for i, s in enumerate(mutations_df.columns.values) if cancer_type == s][0]
    mutations_per_cancer_type[cancer_type] = mutations_df.iloc[:, cancer_type_specific_mut_idx]

    # Split Train/Test
    X_train, X_test, y_train, y_test = train_test_split(scATAC_df,
                                                        mutations_per_cancer_type[cancer_type],
                                                        test_size=0.10,
                                                        random_state=42)

    # Cross Validate
    filename = f"models/{cancer_type}/grid_search_results.pkl"
    os.makedirs(f"models/{cancer_type}", exist_ok=True)
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

    if (not os.path.exists(filename)):
        grid_search_results = grid_search(X_train, y_train, pipe, params, 10)
        pickle.dump(grid_search_results, open(filename, 'wb'))

    grid_search_results = pickle.load(open(filename, 'rb'))

    # Feature Selection
    print(f"Best Score: {grid_search_results.best_score_}\n")
    best_model = grid_search_results.best_estimator_.get_params()['regressor__selected_model']
    n = 20
    top_n_feats = get_top_n_features(best_model, n, scATAC_df.columns.values)
    #print_top_features(top_n_feats)

    backward_eliminate_features(X_train, y_train, best_model, 20, params, 10, cancer_type)

    # Test Set Performance
    test_preds = best_model.predict(X_test)
    test_set_performance = explained_variance_score(y_test, test_preds)
    pickle.dump(test_set_performance, open(f"models/{cancer_type}/test_set_performance.txt", 'wb'))
    print(f"Test set performance: {test_set_performance}")
