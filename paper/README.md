# Reproducing paper figures
The instructions for each figure refer to the lines in the correspondingly named bash script. All data was pre-processed as explained in the README in the homepage of this repo. 

## Figure 1 

### 1a
To build the models for that were used to obtain the COO predictions for Figure 1A, run `../analysis/ML/prep_ML_model_scripts.py` calls below the Figure 1C comment. These will also create bash scripts containing commands for plotting the results (in `../analysis/ML/robustness_scripts`). Each script will be named with a unique ID. Run these from `../analysis/ML` e.g assuming the unique ID is `robustness_c9e9655b-a9b5-4462-a064-db7f69e33ec7`, run 

### 1b
To build the models for that were used to obtain the COO predictions for Figure 1B, run `../analysis/ML/prep_ML_model_scripts.py` calls below the Figure 1B comment (note that you will need to set things up to run jobs in parallel, as explained in the README of the homepage of this repo). To plot the results after building the models, run the call to `../analysis/ML/plot_top_features.R` . This will create a PDF in `figures` called `grid_analysis.pdf`.



```
sh robustness_scripts/robustness_c9e9655b-a9b5-4462-a064-db7f69e33ec7.sh
``` 

from within `../analysis/ML`. 

This will create the figures in 

```
../figures/models/XGB/<CANCER_TYPE>/scATAC_source_<DATASETS>_cell_number_filter_<CELL_NUMBER_FILTER>_annotation_<ANNOTATION>_seed_all_seeds/backwards_elimination_results/
```

in a file called
 
```
<CANCER_TYPE>_feature_importance_with_<TOP_FEATURES_TO_PLOT_FEAT_IMP>_features_top_5_features.pdf
```

### 1c

Instructions analogous to those for 1a.

### 1d
Run the corresponding commands in the Jupyter notebook `../analysis/ML/paper_umaps.ipynb`

### 1e
Run the corresponding commands in the Jupyter notebook `../analysis/ML/fig1e.ipynb`

### 1f
Instructions analogous to those for 1a.

### 1g
Instructions analogous to those for 1a.

### 1h
Run the command under the Figure 1h comment.





## Aside
The next figures require a shared ArchR object. The following steps require quite a bit of memory, so if it doesn't work when you first run it, increase the amount of memory till it works. Also note that for the Shendure dataset, we only need Pancreas and Stomach files. So you can delete the other files or move them to a temporary location. To create this object, we first create Arrow files (can increase number of cores depending on available resources).
```
Rscript ../data/scripts/create_arrow_files_and_tss.R --dataset=Greenleaf_colon --cores=1
Rscript ../data/scripts/create_arrow_files_and_tss.R --dataset=Shendure --cores=1
Rscript ../data/scripts/create_arrow_files_and_tss.R --dataset=Greenleaf_brain --cores=1
Rscript ../data/scripts/create_arrow_files_and_tss.R --dataset=Greenleaf_pbmc_bm --cores=1
```

Once the Arrow files are created, we create a shared ArchR object (again, can increase number of cores depending on available resources):
```
Rscript ../data/scripts/create_ArchR_project.R --cores=1
```

## Figure 2
### 2A
Run the commands below the Figure 2A comment. The first command creates the ArchR object associated with the colon UMAP. It also gets metacells that are needed for the metacell correlation analysis. The next command creates the binned scATAC fragment data per cell, which is needed for obtaining correlations on a per-cell basis (and then using these to obtain metacell correlations). The final command performs the metacell correlation analysis and plots the results in `../figures/mss.pdf`

### 2B
Instructions analogous to those for 1C. 

### 2C
Instructions analogous to those for 2A. 

### 2D
Instructions analogous to those for 1C.

## Figure 3
### 3A
Instructions analogous to those for 1D.

### 3B
Instructions analogous to those for 1C.

### 3C
Instructions analogous to those for 1D. 

### 3D
Instructions analogous to those for 2A. 

### 3E
Instructions analogous to those for 1C. 

### 3F
Instructions analogous to those for 1C. 


## Figure 4
### 4A
Instructions analogous to those for 1C. 

### 4B
Instructions analogous to those for 2A.

### 4C
Instructions analogous to those for 1B. 

### 4D
Instructions analogous to those for 1C. 

### 4E
Start by running `../data/scripts/create_subsampled_data.sh` (uncomment the lines in that script for which you want to create the subsampled data for). Then, instructions are analogous to those for 1C. 

## Supplementary Fig 1
Instructions analogous to those for 1C. Note that the extra figures (top feature appearances and test set performance boxplots) can be found in the folder 

```
../figures/models/XGB/<CANCER_TYPE>/scATAC_source_<DATASETS>_cell_number_filter_<CELL_NUMBER_FILTER>_annotation_<ANNOTATION>_seed_all_seeds/backwards_elimination_results/
```

as

```
<CANCER_TYPE>_top_feature_appearances.pdf
``` 

and 

```
<CANCER_TYPE>_top_feature_test_set_perf_with_<TOP_FEATURES_TO_PLOT>_features.pdf
```

respectively. 

## Supplementary Fig 2
### Sup Fig 2A
Instructions analogous to those for 1C. 

### Sup Fig 2B
For this figure, you need the download the MAF file with the mutation data, which can be found at https://drive.google.com/file/d/1loUifgMF9YkY5_uQ0ARmRZH5hX3po4-E/view?usp=drive_link. Once that's downloaded, put it in `../data/mutation_data` and gunzip the file.
Run the command below the Supplementary Fig2B comment. This will create a figure in `../figures` called `oncoplot.pdf`. 

## Supplementary Fig 3
Instructions analogous to those for 1C. 

## Supplementary Fig 4
Instructions analogous to those for 1C. 

## Supplementary Fig 5
Instructions analogous to those for 1C. 

## Supplementary Fig 6
### Sup Fig 6A
Instructions analogous to those for 1C. 

### Sup Fig 6B
Run the command below the Supplementary Fig6B comment. This will create the marker gene feature plots in `../figures` as `temp.pdf`.  

## Supplementary Fig 7
Instructions analogous to those for 1C. 



