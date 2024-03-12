# Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape

This repository houses the codebase for the paper **Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape**. In this paper, we proposed COCOON (**C**hr**O**matin a**C**cessibility inferred **C**ell Of **O**rigi**N**), which offers a straightforward and cost-effective framework enabling cancer biologists to ascertain the COO at cell-type resolution for a given cancer mutation profile via scATAC sequencing of suspected normal tissues of origin. 
 
In the first section below, we provide an example of how to pre-process data and prepare it as input for COCOON, in addition to then running COCOON. To reproduce paper figures, see section **Reproduce Figures**.

# Example
To run COCOON, we first need to create aggregated, binned scATAC and mutation profiles. We'll start with scATAC, and we'll use the data from [**Single-Cell Multiomic Analysis Identifies Regulatory Programs in Mixed-Phenotype Acute Leukemia**](https://www.nature.com/articles/s41587-019-0332-7) as an example, which comes from PBMC and bonemarrow. After cloning this repository,

`cd data/scripts`

and download the fragment files:

`./get_scATAC_data_from_links.sh ../greenleaf_blood_bone_marrow/greenleaf_blood_bm_ftp_links.txt ../bed_files/greenleaf_pbmc_bm/migrated_to_hg19/ *gz`

This will download the files to the directory `data/bed_files/greenleaf_pbmc_bm/migrated_to_hg19`. 

> **Small aside**: Note that the fragments from this study are already aligned to hg19. However, if your data is not and you would like to run COCOON using the mutation data we used, you will have to align the fragments to hg19, since the mutation data we had available was aligned to hg19. To do this, make sure your fragment files are in the directory `data/bed_files/[DATASET_NAME]/`. Here, DATASET_NAME is whatever name you want to give to your dataset. 

Then run `Rscript data/scripts/migrate_and_save_fragments.R --dataset [DATASET_NAME] --cores [NUM_CORES]`.

We will then rename the files to match the pattern [TISSUE_TYPE]-[SAMPLE_NAME].bed.gz

`cd data/bed_files`

`for f in $(ls *); do echo $(echo $f | sed 's/.*_scATAC_//') | sed 's/.fragments.tsv.gz/.bed.gz/'; done`



