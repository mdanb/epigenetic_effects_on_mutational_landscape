import pandas as pd
import os
import glob
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import numpy as np
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from plotnine import ggplot, geom_point, aes, ggtitle, \
                     labs, theme, scale_color_brewer
from sklearn.preprocessing import StandardScaler

def plot_dendrogram(model, plot_title, linkage, affinity, **kwargs):
    # Create linkage matrix and then plot the dendrogram
    # create the counts of samples under each node
    figure(figsize=(30, 28))
    plt.title(plot_title)

    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)
    plt.xticks(fontsize=13)
    plt.show()
    plt.savefig(f"../figures/cancer_cohort_{plot_title}_linkage_{linkage}_affinity_{affinity}.png")

def save_clustered_datafiles(model, data, cancer_type, distance_threshold, linkage, affinity):
    os.makedirs("../data/processed_data/hierarchically_clustered_mutations/", exist_ok=True)
    threshold_dir = f"../data/processed_data/hierarchically_clustered_mutations/{cancer_type}/linkage_{linkage}_affinity_" \
                    f"{affinity}/dist_thresh_{distance_threshold}"
    os.makedirs(threshold_dir,
                exist_ok=True)
    for idx, label in enumerate(np.unique(model.labels_)):
        patient_cluster = data.loc[model.labels_ == label]
        patient_cluster.to_csv(f"{threshold_dir}/C{idx}.csv")
        patient_cluster.sum(axis=0).to_csv(f"{threshold_dir}/aggregated_C{idx}.csv")

    data.sum(axis=0).to_csv(f"{threshold_dir}/aggregated_all.csv")


def cluster_data_and_display_clustering(distance_threshold, mutations_df, plot_title,
                                        save_clusters=False, affinity="euclidean",
                                        affinity_matrix=None, linkage="ward", no_labels=False):
    model = AgglomerativeClustering(n_clusters=None,
                                    distance_threshold=distance_threshold,
                                    linkage=linkage)
    if (affinity == "euclidean"):
        model = model.fit(mutations_df)
    elif (affinity == "pearson"):
        try:
            model = model.fit(affinity_matrix)
        except:
            print("affinity matrix musn't be None when affinity = 1-pearson")
    plot_dendrogram(model,
                    plot_title,
                    linkage,
                    affinity,
                    truncate_mode="level",
                    p=100,
                    no_labels=no_labels,
                    color_threshold=distance_threshold,
                    above_threshold_color="black",
                    labels = mutations_df.index.to_numpy())
    if (save_clusters):
        save_clustered_datafiles(model, mutations_df, plot_title, distance_threshold,
                                 linkage, affinity)

def get_per_cohort_mutations(cancer_cohort):
    binned_mutation_files = glob.glob(f'../data/processed_data/per_patient_mutations/{cancer_cohort}/binned_mutations*')
    chr_keep = pd.read_csv("../data/processed_data/chr_keep.csv")
    df = pd.DataFrame(columns=chr_keep["chr"].to_numpy())

    for file in binned_mutation_files:
        patient_id = file.split("_")[-1].split(".")[0]
        patient_mutations = pd.read_csv(file)['x']
        df.loc[patient_id] = patient_mutations.values

    agg_mutation_df = pd.read_csv("../data/processed_data/mut_count_data.csv", index_col=0)
    regions_to_consider = agg_mutation_df.loc[~agg_mutation_df.isna().any(axis=1)].index.values
    df = df.loc[:, regions_to_consider]
    return df

def group_donors_by_subtype(mutation_df, annotation_df):
    mutation_df.index = mutation_df.index.to_series().map(
                                                           annotation_df.to_dict()["subtype"]
                                                        )
    return(mutation_df)

def create_cohort_df(grouped_donor_dfs):
    grouped_donor_dfs = pd.concat(grouped_donor_dfs)
    return(grouped_donor_dfs.reset_index(names="subtype"))

def reduce_mutation_dims(df, method, n_components, scale=True):
    if (method == "PCA"):
        model = PCA(n_components=n_components)
    elif (method == "TSNE"):
        model = TSNE(n_components=n_components)
    elif(method == "UMAP"):
        model = umap.UMAP(n_components=n_components, random_state=1)

    # df = df.drop("subtype", axis=1)
    df = np.log(df.div(np.sum(df, axis=1).array, axis=0) + 1)
    if (scale):
        df = StandardScaler().fit_transform(df)
#         df = (df - np.mean(df)) / np.std(df)
#         print(df)
#     print(np.sum(df, axis=0))
#     std_debug = np.std(df, axis=0)
#     mean_debug = np.mean(df, axis=0)

#     print(df.shape)
    df_reduced = model.fit_transform(df)
    return([model, df_reduced])

def plot_reduced_dim_mutations(reduced_dim_data, components_to_plot,
                               mut_counts, plot_title, subtypes):
    df = pd.DataFrame({"dim1": reduced_dim_data[:, 0],
                       "dim2": reduced_dim_data[:, 1],
                       "subtype":subtypes,
                       "mut_counts": mut_counts})
    if (2 in components_to_plot):
        df["dim3"] = reduced_dim_data[:, 2]

    fig, plot = (
        ggplot(df, aes(components_to_plot[0], components_to_plot[1], shape='subtype',
                color=np.log(mut_counts)))
         + geom_point(size=3)
         + ggtitle(plot_title)
         + labs(color="mutation counts")
         + theme(figure_size=(18, 13))
         + scale_color_brewer(type="qual", )
    ).draw(show=False, return_ggplot=True)
    return(fig)

def get_top_variable_bins(df, n):
    stds = np.std(df)
    largest_variance_feat_idxs = np.argsort(stds)[::-1]
    top = df.iloc[:, largest_variance_feat_idxs].iloc[:, 0:n]
    return top

def drop_bottom_count_samples(df, n):
    mut_counts = np.sum(df.drop("subtype", axis=1), axis=1)
    bottom_count_idxs = np.argsort(mut_counts)[0:n]
    top = df.loc[set(df.index) - set(bottom_count_idxs)]
    return top



# Lung
LUAD_US = get_per_cohort_mutations("LUAD-US")
LUSC_US = get_per_cohort_mutations("LUSC-US")
LUAD_US_annotation = pd.DataFrame(["Adeno"] * len(LUAD_US.index), index=LUAD_US.index,
                                 columns = ["subtype"])
LUSC_US_annotation = pd.DataFrame(["Squamous"] * len(LUSC_US.index), index=LUSC_US.index,
                                 columns = ["subtype"])

lung_adeno = group_donors_by_subtype(LUAD_US, LUAD_US_annotation)
lung_squamous = group_donors_by_subtype(LUSC_US, LUSC_US_annotation)
lung = create_cohort_df([lung_adeno, lung_squamous])

bottom_n_remove = 2
lung = drop_bottom_count_samples(lung, bottom_n_remove)
subtype = lung["subtype"]
lung = lung.drop("subtype", axis=1)
n_feats = 2128
lung = get_top_variable_bins(lung, n_feats)
method = "PCA"
model, reduced_dim_data = reduce_mutation_dims(lung, method, 2, True)
fig = plot_reduced_dim_mutations(reduced_dim_data,
                               ["dim1", "dim2"],
                               np.sum(lung, axis=1),
                               "Lung",
                               subtype)
#fig.savefig('fig_log.png', dpi=300)
fig.savefig(f'fig_{method}_top_{n_feats}_features_no_bottom_'
            f'{bottom_n_remove}_samples.png',
            dpi=300)
