# Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape

This repository houses the codebase for the paper **Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape**. In this paper, we proposed COCOON (**C**hr**O**matin a**C**cessibility inferred **C**ell **O**f **O**rigi**N**), which offers a straightforward and cost-effective framework enabling cancer biologists to ascertain the COO at cell-type resolution for a given cancer mutation profile via scATAC sequencing of suspected normal tissues of origin. 
 
In the first section below, we provide an example of how to pre-process data and prepare it as input for COCOON, in addition to then running COCOON. To reproduce paper figures, see scripts in this directory, named by the corresponding figure. 

# Setting up your conda environment
First, we need to up our conda environment. To do this, run

```
sh setup_env.sh
```
Then, activate the environment:

```
conda activate coo
```
and run:
```
conda install postgresql=14.5
conda install -c conda-forge r-gert=2.0.0
conda install -c conda-forge gsl=2.7
```
then:
```
Rscript post_installation.R
```

You should now be ready to proceed.
 
# Example
To run COCOON, we first need to create aggregated, binned scATAC and mutation profiles. All scripts for data processing are in `data/scripts`, and all relative paths below are relative to this directory. 

## scATAC pre-processing
We'll start with scATAC, and we'll use the data from [**Single-Cell Multiomic Analysis Identifies Regulatory Programs in Mixed-Phenotype Acute Leukemia**](https://www.nature.com/articles/s41587-019-0332-7) as an example, which comes from PBMC and bonemarrow. After cloning this repository, download the fragment files:

```
sh get_scATAC_data_from_links.sh ../greenleaf_blood_bone_marrow/greenleaf_blood_bm_ftp_links.txt ../bed_files/Greenleaf_test/migrated_to_hg19/ *gz
```

This will download the files to the directory `data/bed_files/Greenleaf_test/migrated_to_hg19`. 

> **Small aside**: Note that the fragments from this study are already aligned to hg19. However, if your data is not and you would like to run COCOON using the mutation data we used, you will have to align the fragments to hg19, since the mutation data we had available was aligned to hg19. To do this, make sure your fragment files are in the directory `data/bed_files/[DATASET_NAME]/`. Here, `DATASET_NAME` is whatever name you want to give to your dataset. Then run `Rscript migrate_and_save_fragments.R --dataset [DATASET_NAME] --cores [NUM_CORES]`. Note that this script is parallelized, so you can specify `NUM_CORES`, where each core will be migrating one of the fragment files. So if you have 16 fragment files and 4 cores, then specifying `NUM_CORES` to 4 will process 4 files at a time. 

We will then rename the files to match the pattern `[TISSUE_TYPE]-[SAMPLE_NAME].bed.gz`:

```
cd data/bed_files/Greenleaf_test/migrated_to_hg19
```

```
for f in $(ls *); do mv $f $(echo $f | sed 's/.*_scATAC_//' | sed 's/.fragments.tsv.gz/.bed.gz/' | sed 's/_/-/'); done
```

Now, we need a metadata file that contains cell-type annotations. This file must contain the following columns:

- barcode: cell barcode corresponding to the barcodes in the fragment files e.g AACTGGTTCCACTAGA
- sample: sample name from which the cell comes from (this must correspond to `SAMPLE_NAME` in the name of the fragment file, see above). 
- cell_type: cell type annotation
- tissue: tissue name from which the cell comes from (again, must correspond to `TISSUE_TYPE` in the name of the fragment file, see above) 

The metadata file needs to be in `data/metadata` and must be called `[ANNOTATION].csv` where `ANNOTATION` is your desired name for the cell annotations that you used. For this example, we have a file called `test_annotation.csv`. 

We then create the binned scATAC profiles:

```
Rscript create_count_overlaps.R --dataset [DATASET_NAME] --cores [NUM_CORES] --annotation [ANNOTATION] --which_interval_ranges [INTERVAL_RANGES_NAME]
```

This is script is parallelized, and each core will be processing a given fragment file. Note that here, `DATASET_NAME` and `ANNOTATION` must be the ones you chose previously. `INTERVAL_RANGES_NAMES` specifies the GenomicRanges object (bins) that we will align our fragments to. For this example, we have a file called `data/test_ranges.RData`. This corresponds to the bins we used for our paper. So for `INTERVAL_RANGES_NAME`, we would specify `test_ranges`. Note that this file must be an `.RData` file. 

Thus, assuming we have 4 cores, in this example, we would run:

```
Rscript create_count_overlaps.R --dataset Greenleaf_test --cores 4 --annotation test_annotation --which_interval_ranges test_ranges
```

This will create our binned fragment files in `data/processed_data/count_overlap_data/[ANNOTATION]/` with the name `interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_count_overlaps_[TISSUE_TYPE]-[SAMPLE_NAME].rds`. In addition, we will get metadata files in `data/processed_data/cell_counts_per_sample/[ANNOTATION]/` with the names `cell_counts_interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_count_overlaps_[TISSUE_TYPE]-[SAMPLE_NAME].rds`. These files contain data about the number of cells per cell type. 

Finally, we must combine the counts from identical cell types from different samples. We do this by running:

```
Rscript combine_overlaps.R --datasets [DATASET_NAME] --annotation [ANNOTATION] --which_interval_ranges [INTERVAL_RANGES_NAME]
```

In our example:

```
Rscript combine_overlaps.R --datasets Greenleaf_test --annotation test_annotation --which_interval_ranges test_ranges
```

This creates an aggregated scATAC profile, which can be found at `data/processed_data/count_overlap_data/combined_count_overlaps/[ANNOTATION]/interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_combined_count_overlaps.rds`. In the same directory, another file called `interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_combined_count_overlaps_metadata.rds` will also be created, which has information about the number of cells per cell type and tissue type. 

In this example:
`interval_ranges_test_ranges_Greenleaf_test_combined_count_overlaps.rds` and 
`interval_ranges_test_ranges_Greenleaf_test_combined_count_overlaps_metadata.rds`

## Mutation data (SNV) pre-processing
To create aggregated, binned mutation profiles, we first need a [MAF file](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) with the mutation (SNV) data. For this example, we will use Non-Hodgkin lymphoma (Lymph-BNHL). We have the corresponding MAF file in `data/mutation_data/`, called `Lymph-BNHL_SNV_with_SEX.txt`. Your file, if using your own mutation data, should be called `[CANCER_TYPE]_SNV_with_SEX.txt`. Next, we begin by parsing the file:

```
python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types [CANCER_TYPE]
```

We then intersect the mutations with the bins:

```
python3 3_Intersect_paz_cancertypes.py --cancer_types [CANCER_TYPE]
```

Finally, we aggregate the mutations across samples:

```
python3 4_AssembleCout_paz_Cancergroup.py --cancer_types [CANCER_TYPE]
```

This creates a file `data/processed_data/[CANCER_TYPE].txt` that has the aggregated, binned mutation profile.

Finally, we add custom names to the bins, in addition to adding `CANCER_TYPE` as a column name for our data (this step could have been incorporated in the previous scripts):
```
Rscript align_mutations_to_ranges --cancer_types [CANCER_TYPE]
```

This creates a csv file `data/processed_data/[CANCER_TYPE].csv` that is ready to be input into COCOON.

Putting these together for our example:
```
python3 2_Sorting_MutationFileSex_CancerType.py --cancer_types Lymph-BNHL
python3 3_Intersect_paz_cancertypes.py --cancer_types Lymph-BNHL
python3 4_AssembleCout_paz_Cancergroup.py --cancer_types Lymph-BNHL
Rscript align_mutations_to_ranges --cancer_types Lymph-BNHL
```

## Running COCOON
In our paper, for each cancer type, we ran the model 100 times to obtain robust predictions. In practice, this means we needed access to a compute cluster to parallelize the model training process. Below, we present two pipelines for obtaining predictions: an unparallelized approach and a parallelized approach. Of course, to obtain 100 predictions in a reasonable amount of time, particularly if using large feature spaces, you would want to use the parallelized option. In the case of the parallelized option, we assume you have access to a system that uses UGE as the job manager. However, even if this is not the case, it is very straighforward to modify `prep_ML_model_scripts.py`, the script which does the parallelization, to use a different job manager. 

The script `analysis/ML/build_ML_model.py` runs COCOON. Before looking at an example, we outline the important command line options below:

<details>
  <summary><b>Cancer Types [--cancer_types]</b></summary>
  The cancer types to run the model on. Note that this is not parallelized i.e specifying mutliple cancer types will run the model sequentially on each cancer type. Should correspond to <code>CANCER_TYPE</code> above. 
</details>
<details>
  <summary><b>scATAC Datasets [--datasets]</b></summary>
  Names of the scATAC datasets to consider. Should correspond to <code>DATASET_NAME</code> above. 
</details>
<details>
  <summary><b>Cell Type Number Filter [--scATAC_cell_number_filter]</b></summary>
  Filter for the minimum number of cells per cell type. All cells with lower than this minimum will be excluded. 
</details>
<details>
  <summary><b>Annotation [--annotation_dir]</b></summary>
   Name of annotation used. Should correspond to <code>ANNOTATION</code> above. 
</details>
<details>
  <summary><b>Fold For Test Set [--fold_for_test_set]</b></summary>
   Which fold among the 10 contiguous genome regions to use as the test set. This number must be in the range 1-10. 
</details>
<details>
  <summary><b>Seed Range [--seed_range]</b></summary>
   Range of seeds to use when training the model and computing permutation importance for a particular train/test split (i.e for a particular chosen fold for the test set). Format: [LOWEST_SEED]-[HIGHEST_SEED] e.g 3-7, where 3 and 7 are inclusive.
</details>
<details>
  <summary><b>Number Of Optuna Trials Pre-backward Feature Selection [--n_optuna_trials_prebackward_selection]</b></summary>
   Number of Optuna trials to run before running backward feature selection (i.e when using the full feature space)
</details>
<details>
  <summary><b>Number of Optuna Trials During Backward Feature Selection [--n_optuna_trials_backward_selection]</b></summary>
   Number of Optuna trials to run when performing backward feature selection, for each iteration of BFS.
</details>
<details>
  <summary><b>Enable SQLite [--sqlite]</b></summary>
   If specified, this will enable the use of SQLite to store Optuna results. This is the easiest way to store results, as it requires no database configuration on your part. It is only appropriate to set this option when you are not using parallelization. However, SQLite is not well suited for concurrent database access, so it should not be used when training models in parallel. 
</details>
<details>
  <summary><b>Test Set Performance Number Of Features [--test_set_perf_num_features]</b></summary>
   Number of features during BFS at which to save test set performance. This can be a list of numbers e.g 1 2 5 10 or simply "all." It is advisable to just set this option to "all", as computing test set performance requires negligeable compute. However, assuming the models trained during BFS are not deleted, these can always be calculated later (and the script does provide an option to do this by just re-running it, which won't re-train models, but load already trained models and compute the test set performance at the specified number of features). 
</details>
<details>
  <summary><b>Cell Types To Keep [--cell_types_keep]</b></summary>
  Cell types to keep despite not meeting the cell type filter. Format: [NAME_OF_CELL_TYPE]-[DATASET]. It is possible to provide more than 1 e.g C1-D1 C2-D2
</details>
<details>
  <summary><b>Enable Use Of Custom Mutations [--custom_mutations]</b></summary>
  This is one of multiple mutually exclusive options, which specifies where to load mutations from. It should always be set when using your own mutations, and for the example below. Which options to set use when reproducing the results from the paper will be specified as needed.
</details>
<details>
  <summary><b>Tissues To Consider [--tissues_to_consider]</b></summary>
  If you only want to use a particular subset of tissues from your scATAC data, then set this option to the name of these tissues as a list (and exactly as they are called in the metadata file created using <code>combine_overlaps.R</code> under the column "tissue")
</details>
<details>
  <summary><b>Which Interval Ranges [--which_interval_ranges]</b></summary>
  Which interval ranges to use. Should correspond to <code>WHICH_INTERVAL_RANGES</code> used previously. 
</details>

### Unparallelized
We will use the PBMC/bonemarrow + Lymph-BNHL data we generated previously as an example. Run:

```
python3 build_ML_model.py 
--cancer_types Lymph-BNHL
--datasets Greenleaf_test
--scATAC_cell_number_filter 100
--annotation_dir test_annotation
--seed_range 1-1
--fold_for_test_set 1
--n_optuna_trials_prebackward_selection 50
--n_optuna_trials_backward_selection 50
--sqlite
--test_set_perf_num_features all
--custom_mutations
--which_interval_ranges test_ranges
```

Once the model is done training, we can then obtain a plot of the feature importance of different features at a given iteration of BFS using the script `plot_top_features.R`. Before we do this, we outline the important command line options to do this:
<details>
  <summary><b>Cancer Types [--cancer_types]</b></summary>
  The cancer types to plot. Should correspond to cancer type used when training the model. Should be a string, separating tissues by a comma and no space in between e.g "cancer_type1,cancer_type2,..."
</details>
<details>
  <summary><b>Cell Type Number Filter [--cell_number_filter]</b></summary>
  Should correspond to the filter used when training the model.
</details>
<details>
  <summary><b>Datasets [--datasets]</b></summary>
  Should correspond to the datasets used when training the model. Should be a string, separating datasets by a comma and no space in between.
</details>
<details>
  <summary><b>Top features to plot [--top_features_to_plot]</b></summary>
  Number of features at which to plot feature importances (assuming `robustness_analysis` is not enabled, see parallelized option below). Should be a string, separating numbers by a comma and no space in between.
</details>
<details>
  <summary><b>Tissues To Consider [--tissues_to_consider]</b></summary>
  Should correspond to tissues considered when training the model. Should be a string, separating tissues by a comma and no space in between.
</details>
<details>
  <summary><b>Annotation [--annotation]</b></summary>
  Should correspond to annotation used when training the model.
</details>
<details>
  <summary><b>Seed Range [--seed_range]</b></summary>
   Should correspond to seeds used when training the model.
</details>
<details>
  <summary><b>Folds For Test Set [--folds_for_test_set]</b></summary>
   Unlike when building the model, this parameter is "folds" plural i.e it should be a range in the format "A-B" where A, B are integers between 1 and 10 with A <= B
</details>
<details>
  <summary><b>Cell Types To Keep [--cell_types_keep]</b></summary>
  Should correspond to cell types kept when training the model. Should be a string, separating cell types by a comma and no space in between
</details>

For our example, we run:
```
Rscript plot_top_features.R
--cancer_types=Lymph-BNHL
--datasets=Greenleaf_test
--cell_number_filter=100
--annotation=test_annotation
--seed_range=1-1
--top_features_to_plot=10,5,2,1
--folds_for_test_set=1-1
--tissues_to_consider=all
```

This will generate the plot in 
```
figures/models/XGB/Lymph-BNHL/scATAC_source_Greenleaf_test_cell_number_filter_100_annotation_test_annotation_seed_1_fold_for_test_set_1/backwards_elimination_results/
```

in the file called `permutation_importance_bar_plot.png`. The plot is shown below. Note that as we would expect for Lymph-BNHL, the most important feature turns out to be bone marrow B cells (BMMC B), which corresponds to the putative cell of origin. Note that each feature in the plot is appended with the number of features being considered. Since we specified `top_features_to_plot=10,5,2,1`, we see the permutation importance of different features for each of these different number of features left. 

![alt text](https://github.com/mdanb/epigenetic_effects_on_mutational_landscape/blob/main/permutation_importance_bar_plot.png)

### Parallelized
To obtain more confidence in our prediction, we can run the model multiple times, using different random seeds, and different test sets. This is where parallelization comes into play. SQLite, which is what we used in the unparallelized option for storing results from Optuna, is not well suited for concurrent database access. Thus, we used a PostgreSQL database. The first step in thus to set this up. Begin by initializing a new PostgreSQL directory that we'll call `sqldb`:

```
initdb -D sqldb
```
We will change the configurations of the database to make it less restrictive so that our database can accept connections from other hosts (the jobs that will be running the model). To do this, edit the file `analysis/ML/sqldb/postgresql.conf` so that:
- listen_addresses = '*'
- max_connections = 10000

Make sure `listen_addresses` is uncommented. 
Also, in `analysis/ML/sqldb/postgresql.conf`, change the lines under the comments `# IPv4 local connections:` and `# IPv6 local connections:` to:
```
# IPv4 local connections:
host    all             all             0.0.0.0/0            trust
# IPv6 local connections:
host    all             all             ::/0                 trust
```


Next, we need to start the server. We will also get the hostname of the machine that the server is running on; this is helpful for later. Depending on how you do this, the hostname will change every time you start a new server, which you need to do every time you want to run models. Thus, it's best to just create a mini-script that `echo`s the hostname to a file, which will then access from the script that runs the model, and then starts the server:

Create a script that you can run. I'll call it `startpg` (in my case, it's in `~/.local/bin/`, which is one of the directories included in `$PATH`, so running `startpg` will run that script):
  
```
echo $(hostname) > analysis/ML/postgresql_hostname.txt
~/conda/coo/bin/pg_ctl -D analysis/ML/sqldb start
```

You will need to modify the path `~/conda/coo/bin/` depending on where your conda environment is installed.

Now, run:
```
startpg
```

We've now started the server (and we've recorded the hostname of the machine running it). 

We now create a new user for the database called `my_user` (can change to whatever you want). We will be asked for a passcode, and I'll use `password` here:
```
createuser -P my_user
```

We then create a database called `optuna_db` associated with this user (again, we can pick any name):

```
createdb -O my_user optuna_db
```
Now, change line 828 in `ML_utils.py` to reflect the name of the database, as well as the username and password that you chose for your user (note that we're writing the password in the script, which is not good practice, but for our purposes it's fine). In our case:

```
postgresql://my_user:password@{hostname}:5432/optuna_db
```

We're now ready to submit multiple jobs in parallel to train our models. To do this, we will use the script `prep_ML_model_scripts.py`. We first go over the important command line options:
<details>
  <summary><b>Cancer Types [--cancer_types]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>scATAC Datasets [--datasets]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Cell Type Number Filter [--scATAC_cell_number_filter]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Annotation [--annotation_dir]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Test Set Performance Number Of Features [--test_set_perf_num_features]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Number Of Optuna Trials Pre-backward Feature Selection [--n_optuna_trials_prebackward_selection]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Number of Optuna Trials During Backward Feature Selection [--n_optuna_trials_backward_selection]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Cell Types to Keep [--cell_types_keep]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Tissues To Consider [--tissues_to_consider]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Enable Use Of Custom Mutations [--custom_mutations]</b></summary>
  Analogous to option from <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Folds For Test Set [--fold_for_test_set_range]</b></summary>
  Analogous to option from <code>plot_top_features.R</code>, <b>NOT</b> <code>build_ML_model.py</code>
</details>
<details>
  <summary><b>Seeds Interval [--seed_interval]</b></summary>
  Which seeds to use for training. Format: "A-B" where A, B are integers. 
</details>
<details>
  <summary><b>Seed Interval Step [--seed_interval_step]</b></summary>
  Each fold specified by <code>fold_for_test_set_range</code> will be submitted as an individual job. This parameter determines how many seeds per fold to train for each job. So if e.g <code>seed_interval</code> is set to 1-10 and  <code>seed_interval_step</code> is set to 2, then the jobs will train seeds 1-2, 3-4, 5-6, 7-8, and 9-10 respectively, where the seeds within a given range are trained sequentially, not in parallel (so here, if e.g we're considering fold 1, then the first job will first train a model using seed 1, then once that is done, will train the model with seed 2)
</details>
<details>
  <summary><b>Cores [--cores]</b></summary>
  Number of cores to use per job submission 
</details>
<details>
  <summary><b>Memory per Core [--mem_per_core]</b></summary>
  Amount of memory to use per core (in GB)
</details>
<details>
  <summary><b>Time [--time]</b></summary>
  Max time per job. Format: "hours:minutes:seconds"
</details>
<details>
  <summary><b>Submit Jobs [--submit_jobs]</b></summary>
  Submit the jobs (useful to set off for debugging).
</details>

Note that if you are not using UGE, you will need to modify lines 196-205 and 239 to reflect your job manager specifics. 
We start by testing things out before submitting jobs, so we will first omit `--submit_jobs`:

```
python3 prep_ML_model_scripts.py 
--cancer_types Lymph-BNHL 
--scATAC_cell_number_filter 100 
--annotation_dir test_annotation 
--datasets Greenleaf_test 
--seed_interval=1-10 
--n_optuna_trials_prebackward_selection 50 
--n_optuna_trials_backward_selection 50 
--fold_for_test_set_range 1-10 
--test_set_perf_num_features all 
--cores=8 
--seed_interval_step=5 
--custom_mutations 
--mem_per_core 1 
--which_interval_ranges=test_ranges
```

This creates 20 bash scripts in `analysis/ML` which correspond to the jobs that would be submitted if we were to keep `--submit_jobs`. To make sure the jobs will run, we test out one of the scripts:
```
sh cancer_types_Lymph-BNHL_datasets_Greenleaf_test_scATAC_cell_number_filter_100_annotation_dir_test_annotation_tissues_to_consider_all_seed_range_1-5_fold_for_test_set_1.sh
```

If this runs fine, we can stop it and re-run the above command, this time with the `--submit_jobs` option enabled:
```
python3 prep_ML_model_scripts.py 
--cancer_types Lymph-BNHL 
--scATAC_cell_number_filter 100 
--annotation_dir test_annotation 
--datasets Greenleaf_test 
--seed_interval=1-10 
--n_optuna_trials_prebackward_selection 50 
--n_optuna_trials_backward_selection 50 
--fold_for_test_set_range 1-10 
--test_set_perf_num_features all 
--cores=8 
--seed_interval_step=5 
--custom_mutations 
--mem_per_core 1 
--which_interval_ranges=test_ranges
--submit_jobs
```

Once the jobs end, we can then plot results using `plot_top_features.R`. Here, we introduce new command line options that are helpful for the parallelized setting. 

<details>
  <summary><b>Robustness Analysis [--robustness_analysis]</b></summary>
  This specifies that we're operating in the parallelized setting, so that we obtain the correct plots. 
</details>
<details>
  <summary><b>Top Features To Plot [--top_features_to_plot]</b></summary>
  This is not a new option, but it means something different here. When <code>--robustness_analysis</code> is enabled, <code>--top_features_to_plot</code> specifies the number of features remaining at which to plot the R^2 of the model when the top appearing feature across runs is actually the top feature for that particular model e.g if Lung AT2 is the top appearing feature across 100 runs of the model, then specifying <code>--top_features_to_plot=1,2,5,10</code> will plot boxplots of the R^2 when 1, 2, 5, and 10 features remain after backward feature selection, and will only plot those models where the top feature is actually Lung AT2. For an example, see Supplementary Figure 1 of the paper, last column. 
</details>
<details>
  <summary><b>Top Features To Plot For Feature Importance [--top_features_to_plot_feat_imp]</b></summary>
  This specifies the number of features remaining at which to plot feature importances (boxplots) for various features across the runs. Only the top 5 most important features will be plotted, where importance is first measured based on the number of times the feature appears, and ties are broken by median feature importance. For an example, see Figure 1, 1D. 
</details>

We now plot:
```
Rscript plot_top_features.R 
--cancer_types Lymph-BNHL 
--datasets Greenleaf_test 
--cell_number_filter 100 
--annotation test_annotation 
--seed_range 1-10 
--top_features_to_plot_feat_imp 1,2,5,10 
--top_features_to_plot 1,2,5,10 
--folds_for_test_set 1-10 
--robustness_analysis 
--tissues_to_consider all
```


