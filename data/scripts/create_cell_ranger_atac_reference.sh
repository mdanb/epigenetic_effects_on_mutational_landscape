#!/bin/bash
#/seq/regev_genome_portal/SOFTWARE/10X/cellranger-3.1.0/cellranger mkgtf /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/Homo_sapiens.GRCh37.87.gtf /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/Homo_sapiens.GRCh37.87.filtered.gtf \
#                 --attribute=gene_biotype:protein_coding \
#                 --attribute=gene_biotype:lincRNA \
#                 --attribute=gene_biotype:antisense

cd /broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/
/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/software/cellranger-atac-2.1.0/cellranger-atac mkref --config=GRCh37.config 
