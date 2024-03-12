# Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape

This repository houses the codebase for the paper **Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape**. In this paper, we proposed COCOON (**C**hr**O**matin a**C**cessibility inferred **C**ell Of **O**rigi**N**), which offers a straightforward and cost-effective framework enabling cancer biologists to ascertain the COO at cell-type resolution for a given cancer mutation profile via scATAC sequencing of suspected normal tissues of origin. 
 
In the first section below, we provide an example of how to pre-process data and prepare it as input for COCOON, in addition to then running COCOON. To reproduce paper figures, see section **Reproduce Figures**.

# Example
To run COCOON, we first need to create aggregated, binned scATAC and mutation profiles. We'll start with scATAC, and we'll use the data from [**Single-Cell Multiomic Analysis Identifies Regulatory Programs in Mixed-Phenotype Acute Leukemia**](https://www.nature.com/articles/s41587-019-0332-7) as an example, which comes from PBMC and bonemarrow. After cloning this repository,

`cd data/scripts`

and download the fragment files:

`./get_scATAC_data_from_links.sh ../greenleaf_blood_bone_marrow/greenleaf_blood_bm_ftp_links.txt ../bed_files/greenleaf_pbmc_bm/migrated_to_hg19/ *gz`

This will download the files to the directory `data/bed_files/greenleaf_pbmc_bm/migrated_to_hg19`. 

> **Small aside**: Note that the fragments from this study are already aligned to hg19. However, if your data is not and you would like to run COCOON using the mutation data we used, you will have to align the fragments to hg19, since the mutation data we had available was aligned to hg19. To do this, make sure your fragment files are in the directory `data/bed_files/[DATASET_NAME]/`. Here, `DATASET_NAME` is whatever name you want to give to your dataset. Then run `Rscript data/scripts/migrate_and_save_fragments.R --dataset [DATASET_NAME] --cores [NUM_CORES]`. Note that this script is parallelized, so you can specify `NUM_CORES`, where each core will be migrating one of the fragment files. So if you have 16 fragment files and 4 cores, then specifying `NUM_CORES` to 4 will process 4 files at a time. 

We will then rename the files to match the pattern `[TISSUE_TYPE]-[SAMPLE_NAME].bed.gz`:

`cd data/bed_files/greenleaf_pbmc_bm/migrated_to_hg19`

`for f in $(ls *); do echo $(echo $f | sed 's/.*_scATAC_//') | sed 's/.fragments.tsv.gz/.bed.gz/'; done`

Now, we need a metadata file that contains cell-type annotations. This file must contain the following columns:

- barcode: cell barcode corresponding to the barcodes in the fragment files e.g AACTGGTTCCACTAGA
- sample: sample name from which the cell comes from (this must correspond to `SAMPLE_NAME` in the name of the fragment file, see above). 
- cell_type: cell type annotation
- tissue: tissue name from which the cell comes from (again, must correspond to `TISSUE_TYPE` in the name of the fragment file, see above) 

The metadata file needs to be in `data/metadata` and must be called `[ANNOTATION].csv` where `ANNOTATION` is your desired name for the cell annotations that you used. For this example, we have a file called `test_annotation.csv`. 

We then create the binned scATAC profiles:
`Rscript data/scripts/create_count_overlaps.R` --dataset [DATASET_NAME] --cores [NUM_CORES] --annotation [ANNOTATION] --which_interval_ranges [INTERVAL_RANGES_NAMES]


