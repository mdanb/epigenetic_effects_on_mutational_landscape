# Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape

This repository houses the codebase for the paper **Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape**. In this paper, we proposed COCOON (**C**hr**O**matin a**C**cessibility inferred **C**ell Of **O**rigi**N**), which offers a straightforward and cost-effective framework enabling cancer biologists to ascertain the COO at cell-type resolution for a given cancer mutation profile via scATAC sequencing of suspected normal tissues of origin. 
 
In the first section below, we provide an example of how to pre-process data and prepare it as input for COCOON, in addition to then running COCOON. To reproduce paper figures, see scripts in this directory, named by the corresponding figure. 

# Example
To run COCOON, we first need to create aggregated, binned scATAC and mutation profiles. 

## scATAC pre-processing
We'll start with scATAC, and we'll use the data from [**Single-Cell Multiomic Analysis Identifies Regulatory Programs in Mixed-Phenotype Acute Leukemia**](https://www.nature.com/articles/s41587-019-0332-7) as an example, which comes from PBMC and bonemarrow. After cloning this repository,

```
cd data/scripts
```

and download the fragment files:

```
./get_scATAC_data_from_links.sh ../greenleaf_blood_bone_marrow/greenleaf_blood_bm_ftp_links.txt ../bed_files/greenleaf_pbmc_bm/migrated_to_hg19/ *gz
```

This will download the files to the directory `data/bed_files/greenleaf_pbmc_bm/migrated_to_hg19`. 

> **Small aside**: Note that the fragments from this study are already aligned to hg19. However, if your data is not and you would like to run COCOON using the mutation data we used, you will have to align the fragments to hg19, since the mutation data we had available was aligned to hg19. To do this, make sure your fragment files are in the directory `data/bed_files/[DATASET_NAME]/`. Here, `DATASET_NAME` is whatever name you want to give to your dataset. Then run `Rscript data/scripts/migrate_and_save_fragments.R --dataset [DATASET_NAME] --cores [NUM_CORES]`. Note that this script is parallelized, so you can specify `NUM_CORES`, where each core will be migrating one of the fragment files. So if you have 16 fragment files and 4 cores, then specifying `NUM_CORES` to 4 will process 4 files at a time. 

We will then rename the files to match the pattern `[TISSUE_TYPE]-[SAMPLE_NAME].bed.gz`:

```
cd data/bed_files/greenleaf_pbmc_bm/migrated_to_hg19
```

```
for f in $(ls *); do mv $f $(echo $f | sed 's/.*_scATAC_//' | sed 's/.fragments.tsv.gz/.bed.gz/'); done
```

Now, we need a metadata file that contains cell-type annotations. This file must contain the following columns:

- barcode: cell barcode corresponding to the barcodes in the fragment files e.g AACTGGTTCCACTAGA
- sample: sample name from which the cell comes from (this must correspond to `SAMPLE_NAME` in the name of the fragment file, see above). 
- cell_type: cell type annotation
- tissue: tissue name from which the cell comes from (again, must correspond to `TISSUE_TYPE` in the name of the fragment file, see above) 

The metadata file needs to be in `data/metadata` and must be called `[ANNOTATION].csv` where `ANNOTATION` is your desired name for the cell annotations that you used. For this example, we have a file called `test_annotation.csv`. 

We then create the binned scATAC profiles:

```
Rscript data/scripts/create_count_overlaps.R --dataset [DATASET_NAME] --cores [NUM_CORES] --annotation [ANNOTATION] --which_interval_ranges [INTERVAL_RANGES_NAME]
```

This is script is parallelized, and each core will be processing a given fragment file. Note that here, `DATASET_NAME` and `ANNOTATION` must be the ones you chose previously. `INTERVAL_RANGES_NAMES` specifies the GenomicRanges object (bins) that we will align our fragments to. For this example, we have a file called `data/test_ranges.RData`. This corresponds to the bins we used for our paper. So for `INTERVAL_RANGES_NAME`, we would specify `test_ranges`. Note that this file must be an `.RData` file. 

Thus, assuming we have 4 cores, in this example, we would run:

```
Rscript data/scripts/create_count_overlaps.R --dataset Greenleaf_pbmc_bm --cores 4 --annotation test_annotation --which_interval_ranges test_ranges
```

This will create our binned fragment files in `data/processed_data/count_overlap_data/[ANNOTATION]/` with the name `interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_count_overlaps_[TISSUE_TYPE]-[SAMPLE_NAME].rds`. In addition, we will get metadata files in `data/processed_data/cell_counts_per_sample/[ANNOTATION]/` with the names `cell_counts_interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_count_overlaps_[TISSUE_TYPE]-[SAMPLE_NAME].rds`. These files contain data about the number of cells per cell type. 

Finally, we must combine the counts from identical cell types from different samples. We do this by running:

```
Rscript data/scripts/combine_overlaps.R --datasets [DATASET_NAME] --annotation [ANNOTATION] --which_interval_ranges [INTERVAL_RANGES_NAME]
```

This creates an aggregated scATAC profile, which can be found at `data/processed_data/count_overlap_data/combined_count_overlaps/[ANNOTATION]/interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_combined_count_overlaps.rds`. In the same directory, another file called `interval_ranges_[INTERVAL_RANGES_NAME]_[DATASET_NAME]_combined_count_overlaps_metadata.rds` will also be created, which has information about the number of cells per cell type and tissue type. 

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

This creates a file `data/processed_data/[CANCER_TYPE].txt`, or `data/processed_data/Lymph-BNHL.txt` in this case, that has the aggregated, binned mutation profile.

Finally, we add custom names to the bins, in addition to adding `CANCER_TYPE` as a column name for our data (this step could have been incorporated in the previous scripts):
`Rscript align_mutations_to_ranges --cancer_types [CANCER_TYPE]`

This creates a csv file `data/processed_data/[CANCER_TYPE].csv` that is ready to be input into COCOON.

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
  <summary><b>Cell Type Filter [--scATAC_cell_number_filter]</b></summary>
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

### Unparallelized

We will use the PBMC/bonemarrow + Lymph-BNHL data we generated previously as an example. Run:

```
python3 analysis/ML/build_ML_model.py 
--cancer_types Lymph-BNHL
--datasets Greenleaf_pbmc_bm

```

### Parallelized
The con





