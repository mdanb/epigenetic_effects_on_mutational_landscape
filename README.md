# Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape

This repository houses the codebase for the paper **Identifying the origin of cancer at cell-type resolution by modeling the relationship between scATAC genome accessibility and the tumor mutational landscape**. In this paper, we proposed COCOON (**C**hr**O**matin a**C**cessibility inferred **C**ell Of **O**rigi**N**), which offers a straightforward and cost-effective framework enabling cancer biologists to ascertain the COO at cell-type resolution for a given cancer mutation profile via scATAC sequencing of suspected normal tissues of origin. 
 
In the first section below, we provide an example of how to pre-process data and prepare it as input for COCOON, in addition to then running COCOON. To reproduce paper figures, see section **Reproduce Figures**.

# Example
To run COCOON, we first need to create aggregated, binned scATAC and mutation profiles. We'll start with scATAC, and we'll use the data from [**Single-Cell Multiomic Analysis Identifies Regulatory Programs in Mixed-Phenotype Acute Leukemia**](https://www.nature.com/articles/s41587-019-0332-7) as an example, which comes from PBMC and bonemarrow. After cloning this repository,

`cd data/scripts`


