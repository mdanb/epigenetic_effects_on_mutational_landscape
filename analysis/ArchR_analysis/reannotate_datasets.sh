#!/bin/bash


#Rscript reannotate_datasets.R --cores=8 --dataset=Bingren --metadata_for_celltype_fn=GSE184462_metadata.tsv --sep_for_metadata=\t --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=frontal_cortex --nfrags_filter=1000 --tss_filter=4 --cell_types=all --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn="bingren_remove_same_celltype_indexing.csv" --sep_for_metadata="," --cell_type_col_in_metadata="cell.type" --cell_name_col_in_metadata="cellID" --column_to_color_by="cell.type" --tissue="frontal_cortex" --nfrags_filter=1000 --tss_filter=4 --cell_types=all --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Shendure" --metadata_for_celltype_fn="GSE149683_File_S2.Metadata_of_high_quality_cells.txt" --sep_for_metadata="$(echo -e '\t')" --cell_type_col_in_metadata=cell_type --cell_name_col_in_metadata=cell --column_to_color_by=cell_type --tissue=cerebrum --nfrags_percentile=0.2 --tss_percentile=0.2 --cell_types=all

#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn="bingren_remove_same_celltype_indexing.csv" --sep_for_metadata="," --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=frontal_cortex --cell_types=all --marker_genes=OLIG1,OLIG2,AQP4,RBFOX3,EOMES --min_cells_per_cell_type=100





#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn="combined_distal_proximal.csv" --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --cell_name_col_in_metadata=X --column_to_color_by=NULL --tissue=all --nfrags_percentile=0.2 --tss_percentile=0.2 --cell_types="all" --cluster_res=1.5
#Rscript reannotate_datasets.R --cores=8 --dataset="Bingren" --metadata_for_celltype_fn=GSE184462_metadata.tsv --sep_for_metadata="$(echo -e '\t')" --cell_type_col_in_metadata=cell.type --cell_name_col_in_metadata=cellID --column_to_color_by=cell.type --tissue=stomach --cell_types="all" --min_cells_per_cell_type=100 --marker_genes="TFF1,MUC2,TFF3,ATP4A,MUC5B,CLCA1,KLF4,MUC6,FUT2,REG4,AGR2,SPDEF" --min_cells_per_cell_type=100

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --cluster --cluster_res=1 --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --plot_cell_types --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --plot_cell_types --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --cluster --cluster_res=0.3 --plot_cell_types

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --nfrags_percentile=0.2 --tss_percentile=0.2 --plot_cell_types

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --nfrags_percentile=0.15 --tss_percentile=0.15 --plot_cell_types
#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --plot_doublets

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_per_cell_type --plot_doublets

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --plot_doublets

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --cluster --cluster_res=1.5 --filter_doublets --filter_per_cell_type

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --cluster --cluster_res=0.5 --filter_doublets --filter_per_cell_type

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_doublets --filter_per_cell_type --marker_genes="MMP10,KRT5,DLK2,IL33,TSLP,TP63,MKI67,TOP2A,SCGB1A1,SCGB3A1,FOXJ1,FOXN4,TMEM190,GSTA1,GSTA2,ASCL1,CHGA,FOXI1,POU2F3,ASCL2,ENO1,ANXA2,KRT17,AGER,EMP2,SFTPA1,KRT5,KRT17,KRT15,FOXJ1,CAPS,BPIFB1,PIGR,MUC5B,ASCL1,FOXI1,ITLN1,MYH11,COL1A1,COL1A2,ASPN,EPCAM,CALB2,MSLN,MUC1,KRT6A,PECAM1,PTPRC,CD163,LYZ,FCER1G,C1QA,C1QB,C1QC,APOC1,APOE,CD79B,CD3E,CD3D,IL32,CD2,CXCR4,NKG7,GATA4,GATA6,WT1"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_doublets --filter_per_cell_type --marker_genes="CD3D,CD3E,CD3G,CD247,CD4,CD8A,CD8B,IL2RA,FOXP3,IL7R,PTPRC,MS4A1,CD19,NCAM1,ITGAM,CD14,FCGR3A,FCGR3B,FCGR2A,FCGR2B"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=IC --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_doublets --filter_per_cell_type --marker_genes="MMP10,KRT5,DLK2,IL33,TSLP,TP63,MKI67,TOP2A,SCGB1A1,SCGB3A1,FOXJ1,FOXN4,TMEM190,GSTA1,GSTA2,ASCL1,CHGA,FOXI1,POU2F3,ASCL2,ENO1,ANXA2,KRT17,AGER,EMP2,SFTPA1,KRT5,KRT17,KRT15,FOXJ1,CAPS,BPIFB1,PIGR,MUC5B,ASCL1,FOXI1,ITLN1,MYH11,COL1A1,COL1A2,ASPN,EPCAM,CALB2,MSLN,MUC1,KRT6A,PECAM1,PTPRC,CD163,LYZ,FCER1G,C1QA,C1QB,C1QC,APOC1,APOE,CD79B,CD3E,CD3D,IL32,CD2,CXCR4,NKG7,GATA4,GATA6,WT1"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=all --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_per_cell_type --marker_genes="SCGB3A1,SCGB3A2,TP63,KRT5"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=all --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --marker_genes="CD79A,CD79B"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=all --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --de_novo_marker_discovery --cluster_res=0.5 --filter_doublets

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=all --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --marker_genes="CD14,FCGR3A,FCGR3B,FCGR1A,CSF1R,CCR2,CD68,CD163,HLA-DRB1,HLA-DRB5,ITGAM,ITGAX,MRC1,CD1C,THBD,CLEC4C,CD80,CD86"

# Neuronal: POU3F2, NEUROD1, INSM1, NEUROG1, OLIG2, RFX2, RFX5, POU2F1
# Myeloid: C1QA, C1QC, CD68, CD163, THBS1, FN1, C5AR2, S100A8, S100A9, S100A12, FCN1, CD14, VCAN
Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --filter_doublets --filter_per_cell_type --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=all --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --marker_genes="CEBPE,GFI1,IRF1,GNL2,ELANE,PPARG,BHLHE41,EGR2,FABP4,HBA1,MARCO,MME,F13A1,SLC40A1,SEPP,STAB1,FOLR2,CD1C,CD1A,FCER1A,HLA-DQA1,CLEC10A,PKIB,CLEC9A,GCSAM,BATF3,WDFY4,LILRA4,GZMB,CD123,LAMP3,CCR7,FSCN1,CCL22,MARCKSL1,EBI3,IDO1,S100A8,S100A9,S100A12,FCN1,CD14,VCAN,PTX3,LRG1,LYPD2,LST1,LILRB2,FCGR3A,NKG7,NCR1,SPON2,KLRD1,GNLY,KLRC1,KLRF1,FGFBP2,KIT,TPSAB1,CPA3,CTSG,BACH2,ITGA2B,GP1BA,VWF,FLI1,MS4A1,SDC1,KLK1,EBF1,CD79A,CD79B,PAX5,VPREB3,POU2AF1,MZB1,XBP1,CD3D,TCRA,LEF1,TCF7,TCF3,IL7R,CD8A,CD4,BCL11B,FOXP3,IL2RA,TNFRSF4,TNFRSF18,CTLA4,EPCAM,CLU,CLDN18,HNF1B,PECAM1,KDR,PLVAP,MCAM,COL1A1,COL1A2,SPARCL1,PDGFRA,NEUROG1,OLIG2,RFX2,NOTO,POU2F1,RFX5,C1QA,C1QC,CD68,CD163,THBS1,FN1,C5AR2"

#Rscript reannotate_datasets.R --cores=8 --dataset="Tsankov" --metadata_for_celltype_fn=combined_distal_proximal.csv --sep_for_metadata="," --cell_type_col_in_metadata=celltypes --tissue=RPL --cell_types="all" --nfrags_filter=1 --tss_filter=0 --min_cells_per_cell_type=1 --filter_doublets --filter_per_cell_type --marker_genes="ACTA2,COL1A1,COL1A2,VIM,CD44,NT5E,THY1,ENG"
