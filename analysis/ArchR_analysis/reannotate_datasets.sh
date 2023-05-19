#!/bin/bash


#Rscript reannotate_datasets.R --cores=8 --dataset=Bingren --metadata_for_celltype_fn=GSE184462_metadata.tsv --sep_for_metadata=\t --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=frontal_cortex --nfrags_filter=1000 --tss_filter=4 --cell_types=all --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn="bingren_remove_same_celltype_indexing.csv" --sep_for_metadata="," --cell_type_col_in_metadata="cell.type" --cell_name_col_in_metadata="cellID" --column_to_color_by="cell.type" --tissue="frontal_cortex" --nfrags_filter=1000 --tss_filter=4 --cell_types=all --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Shendure" --metadata_for_celltype_fn="GSE149683_File_S2.Metadata_of_high_quality_cells.txt" --sep_for_metadata="$(echo -e '\t')" --cell_type_col_in_metadata=cell_type --cell_name_col_in_metadata=cell --column_to_color_by=cell_type --tissue=cerebrum --nfrags_percentile=0.2 --tss_percentile=0.2 --cell_types=all

#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn="bingren_remove_same_celltype_indexing.csv" --sep_for_metadata="," --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=frontal_cortex --cell_types=all --marker_genes=OLIG1,OLIG2,AQP4,RBFOX3,EOMES --min_cells_per_cell_type=100





#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn="combined_distal_proximal.csv" --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --cell_name_col_in_metadata=X --column_to_color_by=NULL --tissue=all --nfrags_percentile=0.2 --tss_percentile=0.2 --cell_types="all" --cluster_res=1.5
#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn=GSE184462_metadata.tsv --sep_for_metadata="$(echo -e '\t')" --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=stomach --cell_types="all" --min_cells_per_cell_type=100 --marker_genes="TFF1,MUC2,TFF3,ATP4A,MUC5B,CLCA1,KLF4,MUC6,FUT2,REG4,AGR2,SPDEF" --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --cluster --cluster_res=1 --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --plot_cell_types --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1

Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --plot_cell_types --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1

#Rscript reannotate_datasets.R --dataset=Tsankov --cores=8 --cluster --cluster_res=1.2 --marker_genes=WT1,ITLN1,COL1A1,PECAM1,LYZ,CD3D,MSLN,KRT18,KRT5,VIM
