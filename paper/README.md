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
Run the corresponding commands in the Jupyter notebook `../analysis/paper_umaps.ipynb`

### 1e
Run the corresponding commands in the Jupyter notebook `../analysis/fig1e.ipynb`

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
### 2a
Run the commands below the Figure 2A comment. The first command creates the ArchR object associated with the colon UMAP. It also gets metacells that are needed for the metacell correlation analysis. The next command creates the binned scATAC fragment data per cell, which is needed for obtaining correlations on a per-cell basis (and then using these to obtain metacell correlations). The final command performs the metacell correlation analysis and plots the results in `../figures/mss.pdf`

### 2b
Instructions analogous to those for 1a. 

### 2c
Instructions analogous to those for 2c. 

### 2d
Instructions analogous to those for 1a.

## Figure 3
### 3a
Instructions analogous to those for 1a.

### 3b
Instructions analogous to those for 1d.

### 3C
Instructions analogous to those for 1d. 

### 3d
Instructions analogous to those for 2a. 

### 3e
Instructions analogous to those for 1a. 

### 3g
Instructions analogous to those for 1a. 


## Figure 4
### 4a
Run the corresponding commands in the path `../metaplasia_analysis`

### 4b
Run the corresponding commands in the path `../metaplasia_analysis`

### 4c
Run the corresponding commands in the path `../metaplasia_analysis`

### 4d
Run the corresponding commands in the path `../metaplasia_analysis`

### 4e
Instructions analogous to those for 1a. 

## Figure 5

### fig 5a
Instructions analogous to those for 1a.

### fig 5b
Instructions analogous to those for 2a. 

### fig 5c
Instructions analogous to those for 1c. 


### fig 5d
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
### Sup Fig 2a
Instructions analogous to those for 1c. 

### Sup Fig 2b
For this figure, you need the download the MAF file with the mutation data, which can be found at https://drive.google.com/file/d/1loUifgMF9YkY5_uQ0ARmRZH5hX3po4-E/view?usp=drive_link. Once that's downloaded, put it in `../data/mutation_data` and gunzip the file.
Run the command below the Supplementary Fig2b comment. This will create a figure in `../figures` called `oncoplot.pdf`. 

## Supplementary Fig 3
Instructions analogous to those for 1c. 

## Supplementary Fig 4
Instructions analogous to those for 1c. 

## Supplementary Fig 5
Instructions analogous to those for 1c. 

## Supplementary Fig 6
### Sup Fig 6a
Run the corresponding commands in the path `../metaplasia_analysis`

### Sup Fig 6b
Run the corresponding commands in the path `../metaplasia_analysis`

### Sup Fig 6c
Instructions analogous to those for 1c. 

### Sup Fig 6d
Instructions analogous to those for 1c. 

## Supplementary Fig 7
Instructions analogous to those for 1c. 

## Supplementary Fig 8
### Sup Fig 8a
Instructions analogous to those for 1c. 

### Sup Fig 8b
Run the command below the Supplementary Fig6B comment. This will create the marker gene feature plots in `../figures` as `temp.pdf`.  

## Supplementary Fig 9
Instructions analogous to those for 1c. 

## Supplementary Fig 10
Instructions analogous to those for 1c. 

## Supplementary Fig 11
Run the corresponding commands in the Jupyter notebook `../analysis/sup11.ipynb`




