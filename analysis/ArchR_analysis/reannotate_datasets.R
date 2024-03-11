library(ArchR)
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(dplyr)
library(ComplexHeatmap)
library(parallel)
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/knnGen.R')
source("../../utils.R")

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character", default="all"),
  make_option("--metadata_for_celltype_fn", type="character"),
  make_option("--sep_for_metadata", type="character", default=","),
  make_option("--cell_type_col_in_metadata", type="character"),
  # make_option("--cell_name_col_in_metadata", type="character"),
  make_option("--cluster", action="store_true", default=F),
  make_option("--cluster_res", type="double"),
  make_option("--plot_custom_column", action="store_true", default=F),
  make_option("--color_embedding_by", type="character", default="cell_type"),
  make_option("--plus_to_add_to_metadata", type="character", default=NULL),
  make_option("--plus_filters", type="character", default=NULL),
  make_option("--tissue", type="character", default="all"),
  make_option("--nfrags_filter", type="integer", default=1),
  make_option("--tss_filter", type="integer", default=0),
  make_option("--tss_percentile", type="double", default=NULL),
  make_option("--nfrags_percentile", type="double", default=NULL),
  make_option("--filter_per_cell_type", action="store_true", default=FALSE),
  make_option("--cell_types", type="character", default="all"),
  make_option("--marker_genes", type="character", default=NULL),
  make_option("--min_cells_per_cell_type", type="integer", default=1),
  make_option("--plot_doublet_scores", action="store_true", default=FALSE),
  make_option("--filter_doublets", action="store_true", default=FALSE),
  make_option("--save_clusters", action="store_true", default=FALSE),
  make_option("--reannotate", action="store_true", default=FALSE),
  make_option("--de_novo_marker_discovery", action="store_true", default=FALSE),
  # make_option("--plot_correlation_with_cancer", action="store_true", default=FALSE),
  # make_option("--cancer_for_correlation_plot", type="character"),
  make_option("--harmonize", action="store_true", default=FALSE),
  make_option("--doublet_filter", type="double", default=FALSE),
  make_option("--get_metacells", action="store_true", default=FALSE)
)

# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--cores=8",
#                      "--dataset=Tsankov",
#                      "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                      "--sep_for_metadata=,",
#                      "--cell_type_col_in_metadata=celltypes",
#                      "--cluster",
#                      "--cluster_res=0.5",
#                      "--tissue=RPL",
#                      "--nfrags_filter=1",
#                      "--tss_filter=0",
#                      "--cell_types=all",
#                      "--min_cells_per_cell_type=1",
#                      "--filter_per_cell_type",
#                      "--filter_doublets", 
#                      "--marker_genes=MMP10,KRT5,DLK2,IL33,TSLP,TP63,MKI67,TOP2A,SCGB1A1,SCGB3A1,FOXJ1,FOXN4,TMEM190,GSTA1,GSTA2,ASCL1,CHGA,FOXI1,POU2F3,ASCL2,ENO1,ANXA2,KRT17,AGER,EMP2,SFTPA1,KRT5,KRT17,KRT15,FOXJ1,CAPS,BPIFB1,PIGR,MUC5B,ASCL1,FOXI1,ITLN1,MYH11,COL1A1,COL1A2,ASPN,EPCAM,CALB2,MSLN,MUC1,KRT6A,PECAM1,PTPRC,CD163,LYZ,FCER1G,C1QA,C1QB,C1QC,APOC1,APOE,CD79B,CD3E,CD3D,IL32,CD2,CXCR4,NKG7,GATA4,GATA6,WT1")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cluster",
#                       "--cluster_res=0.5",
#                       "--tissue=IC",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--marker_genes=MMP10,KRT5,DLK2,IL33,TSLP,TP63,MKI67,TOP2A,SCGB1A1,SCGB3A1,FOXJ1,FOXN4,TMEM190,GSTA1,GSTA2,ASCL1,CHGA,FOXI1,POU2F3,ASCL2,ENO1,ANXA2,KRT17,AGER,EMP2,SFTPA1,KRT5,KRT17,KRT15,FOXJ1,CAPS,BPIFB1,PIGR,MUC5B,ASCL1,FOXI1,ITLN1,MYH11,COL1A1,COL1A2,ASPN,EPCAM,CALB2,MSLN,MUC1,KRT6A,PECAM1,PTPRC,CD163,LYZ,FCER1G,C1QA,C1QB,C1QC,APOC1,APOE,CD79B,CD3E,CD3D,IL32,CD2,CXCR4,NKG7,GATA4,GATA6,WT1")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cluster",
#                       "--cluster_res=0.5",
#                       "--tissue=IC",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--marker_genes=ACTA2,COL1A1,COL1A2,VIM,CD44,CD73,CD90,CD105")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cluster",
#                       "--cluster_res=0.5",
#                       "--tissue=IC",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--marker_genes=CD3D,CD3E,CD3G,CD247,CD4,CD8A,CD8B,IL2RA,FOXP3,IL7R,PTPRC,MS4A1,CD19,NCAM1,ITGAM,CD14,FCGR3A,FCGR3B,FCGR2A,FCGR2B")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=tsankov_refined_annotation.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--filter_doublets",
#                       "--marker_genes=KRT15,KRT17,KRT5,S100A2,EPCAM,KRT4,KRT13,TP63,SOX2,HES2,FOXA1,SOX4,NKX2-1,MUC5B,SCGB1A1,SCGB3A1,SCGB3A2")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=tsankov_refined_annotation.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=Basal",
#                       "--min_cells_per_cell_type=1",
#                       "--filter_doublets",
#                       "--marker_genes=KRT15,KRT17,KRT5,S100A2,EPCAM,KRT4,KRT13,TP63,SOX2,HES2,FOXA1,SOX4,NKX2-1,MUC5B,SCGB1A1,SCGB3A1,SCGB3A2")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--marker_genes=CD3D,CD3E,CD3G,CD247,CD4,CD8A,CD8B,IL2RA,FOXP3,IL7R,PTPRC,MS4A1,CD19,NCAM1,ITGAM,CD14,FCGR3A,FCGR3B,FCGR2A,FCGR2B")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--marker_genes=KIT,FCER1A,TPSAB1,TPSB2,CPA3,CMA1,HNMT,HRH1")
# )
# 

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--de_novo_marker_discovery",
#                       "--cluster_res=0.6",
#                       "--filter_doublets")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--de_novo_marker_discovery",
#                       "--cluster_res=0.6",
#                       "--filter_doublets")
# )

args = parse_args(OptionParser(option_list=option_list), args=
                    c("--cores=8",
                      "--dataset=Tsankov",
                      "--metadata_for_celltype_fn=tsankov_default_annotation.csv",
                      "--sep_for_metadata=,",
                      "--cell_type_col_in_metadata=celltypes",
                      "--tissue=all",
                      "--nfrags_filter=1",
                      "--tss_filter=0",
                      "--cell_types=all",
                      "--min_cells_per_cell_type=1",
                      "--de_novo_marker_discovery",
                      "--cluster_res=0.6",
                      "--filter_doublets",
                      "--filter_per_cell_type")
)


# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_colon",
#                       "--metadata_for_celltype_fn=greenleaf_colon_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=general_cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=epithelial",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--plus_to_add_to_metadata=GrossPathology,CellType",
#                       "--color_embedding_by=GrossPathology",
#                       "--marker_genes=MSLN,AQP5,TACSTD2,FSCN1,TFF2,ANXA1,ANXA10,REG4,MUC17,S100P,GSDMB,GSDMD,IL18,RELB,MDK,AHR,PDX1"
# ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_colon",
#                       "--metadata_for_celltype_fn=greenleaf_colon_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=general_cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=10000",
#                       "--tss_filter=0",
#                       "--cell_types=epithelial",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--plus_to_add_to_metadata=GrossPathology,CellType",
#                       "--plus_filters=Normal|Unaffected",
#                       "--color_embedding_by=GrossPathology",
#                       # "--marker_genes=CLDN2,CD44,AXIN2,RNF43,TGFBI,EPHB2,TEAD2,CDX2,LGR5,OLFM4,ASCL2",
#                       "--harmonize"
#                     ))
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Yang_kidney",
#                       "--metadata_for_celltype_fn=41467_2021_27660_MOESM4_ESM.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltype",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--color_embedding_by=batch",
#                       "--filter_doublets",
#                       "--doublet_filter=5"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_brain",
#                       "--metadata_for_celltype_fn=GSE162170_atac_cell_metadata.txt.gz",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--get_metacells"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_colon",
#                       "--metadata_for_celltype_fn=greenleaf_colon_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=general_cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=epithelial",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--plus_to_add_to_metadata=GrossPathology,CellType",
#                       "--color_embedding_by=GrossPathology",
#                       "--marker_genes=MDK,ELF3,MSLN,RAB15,CXCL16,ADAM9,HES4,HES1,AQP5,ABHD4,AHNAK,AK1,AKR1B10,ANXA1,ANXA3,BMP8B,BOK,CD55,CLIC3,CRIP2,EPS8L1,DAPK1,DCXR,ECM1,FOSL1,GJB3,GSN,HSPB1,HYAL1,IL1RN,ITGB4,KIFC3,LMNA,PHLDA2,PHLDA3,PDLIM7,P2RY2,PDZK1IP1,PLAUR,PRSS22,CAVIN3,PLCD3,PSCA,RHOD,S100A11,S100A14,S100A16,S100A4,SERPINB5,SLC45A3,TACSTD2,TIMP2,TTC9,VAMP5,VNN1,VSIG1,WWC2"
#                     ))
# general_cell_type

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Rawlins_fetal_lung",
#                       "--metadata_for_celltype_fn=rawlins_fetal_lung_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--color_embedding_by=cell_type"
#                     ))

args = parse_args(OptionParser(option_list=option_list), args=
                    c("--cores=4",
                      "--dataset=Bingren",
                      "--metadata_for_celltype_fn=GSE184462_metadata.tsv",
                      "--sep_for_metadata=\t",
                      "--cell_type_col_in_metadata=celltype",
                      "--tissue=all",
                      "--nfrags_filter=1",
                      "--tss_filter=0",
                      "--min_cells_per_cell_type=1",
                      "--cluster_res=0.6",
                      "--filter_per_cell_type",
                      "--plot_custom_column",
                      "--color_embedding_by=batch"
                    ))
# 
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Shendure",
#                       "--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       "--plot_custom_column",
#                       "--color_embedding_by=batch",
#                       "--marker_genes=ACTA2,TAGLN,MYH11"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_pbmc_bm",
#                       "--metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=10000",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=100",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       # "--plot_custom_column"
#                       "--color_embedding_by=cell_type"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_pbmc_bm",
#                       "--metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=5000",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=100",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       # "--plot_custom_column"
#                       "--color_embedding_by=cell_type"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_pbmc_bm",
#                       "--metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=7500",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=100",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       # "--plot_custom_column"
#                       "--color_embedding_by=cell_type"
#                     ))

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_pbmc_bm",
#                       "--metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=100",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       # "--plot_custom_column"
#                       "--color_embedding_by=cell_type"
#                     ))
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Greenleaf_pbmc_bm",
#                       "--metadata_for_celltype_fn=intermediate_blood_bm_annotation_metadata.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--min_cells_per_cell_type=1",
#                       "--cluster_res=0.6",
#                       "--filter_per_cell_type",
#                       # "--plot_custom_column"
#                       "--color_embedding_by="
#                     ))

# plus means other stuff as well
add_cell_types_plus_to_cell_col_data <- function(cell_col_data, metadata,
                                                 cell_type_col_in_orig_metadata, 
                                                 dataset, plus_to_add_to_metadata) {
  if (dataset == "Shendure") {
    to_match = paste(metadata[["sample_name"]], metadata[["cell"]], sep="#")
    rownames_archr = rownames(cell_col_data)
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    to_match = paste(metadata[["sample"]], metadata[["barcode"]], sep="#")
    rownames_archr = rownames(cell_col_data)
  }
  else if (dataset == "Greenleaf_brain") {
    to_match = paste(metadata[["Sample.ID"]], metadata[["Cell.Barcode"]], sep="#")
    rownames_archr = rownames(cell_col_data)
  }
  else if (dataset == "Yang_kidney") {
    metadata[metadata[["batch"]] == 1, "sample"] = "SRR13679156"
    metadata[metadata[["batch"]] == 2, "sample"] = "SRR13679157"
    to_match = paste(metadata[["sample"]], metadata[["barcode"]], sep="#")
    to_match = gsub("-2", "-1", to_match)
    rownames_archr = rownames(cell_col_data)
  }
  else if (dataset == "Bingren") {
    metadata = metadata[metadata["Life.stage"] == "Adult", ]
    metadata["cellID"] = gsub("\\+", "#", metadata[["cellID"]])
    # idx = str_locate(metadata[["cellID"]], "#")[, 1] 
    # to_match = unname(mapply(function(s, i, ins) {
    #   paste0(substring(s, 1, i-1), ins, substring(s, i))
    # }, s = metadata[["cellID"]], i = idx-1, ins = "rep"))
    to_match = metadata[["cellID"]]
    rownames_archr = rownames(cell_col_data)
    # rownames_archr = unlist(lapply(lapply(strsplit(rownames(cell_col_data), split="_"), 
    #                      "[", 2:4), paste, collapse="_"))
  }
  else if (dataset == "Tsankov") {
    rownames_archr = rownames(cell_col_data)
    to_match = metadata[["X"]]
  }
  else if (dataset == "Greenleaf_colon") {
    # rownames_archr = unlist(lapply(rownames(cell_col_data), 
    #                               function(x) substr(x, 1, nchar(x) - 2)))
    rownames_archr = strsplit(rownames(cell_col_data), split = "_")
    rownames_archr = lapply(rownames_archr, function(x) x[2:length(x)])
    
    exceptions = c("A001-C-104-D_20200214", "A001-C-104-D_20200811",
                   "A001-C-124-D_20200214", "A001-C-124-D_20200702",
                   "A001-C-023-D_20200214", "A001-C-023-D_20200715",
                   "A002-C-010-D_20200310", "A002-C-010-D_20200702",
                   "B004-A-004-D_20200817")
    exceptions = paste(exceptions, collapse = "|")
    exceptions = grep(exceptions, rownames(cell_col_data))
    rownames_archr[exceptions] = lapply(rownames_archr[exceptions],
                                       paste,
                                       collapse="_")
    rownames_archr[-exceptions] = 
      unlist(lapply(rownames_archr[-exceptions], function(x) 
             paste(x[1], unlist(strsplit(x[2], "#"))[2], sep="#")))
    to_match = metadata[["Cell"]]
    rownames_archr = unlist(rownames_archr)
  } else if (dataset == "Rawlins_fetal_lung") {
      WSSS_F_idx = grep("WSSS_F", rownames(cell_col_data))
      WSSS_F_names = lapply(strsplit(rownames(cell_col_data)[WSSS_F_idx], 
                                     split="_"), "[", 4:6)
      WSSS_F_names = unlist(lapply(WSSS_F_names, paste, collapse="_"))
      other = strsplit(rownames(cell_col_data)[-WSSS_F_idx], split="_")
      length_other = lapply(other, length)
      l_6_idx = length_other == 6
      l_6_names = unlist(lapply(other[l_6_idx], "[", 3))
      l_8_idx = length_other == 8
      l_8_names = unlist(lapply(other[l_8_idx], "[", 4))
      other[l_6_idx] = l_6_names
      other[l_8_idx] = l_8_names
      other = unlist(other)
      cell_ids = unlist(lapply(strsplit(rownames(cell_col_data), "#"), "[", 2))
      rownames_archr = unlist(lapply(strsplit(rownames(cell_col_data), "#"), "[", 1))
      rownames_archr[WSSS_F_idx] = WSSS_F_names
      rownames_archr[-WSSS_F_idx] = other
      rownames_archr = paste(rownames_archr, cell_ids, sep="#")
      to_match = metadata[["X"]]
  }
  idx_cell_id = match(rownames_archr, to_match)
  cell_types = metadata[[cell_type_col_in_orig_metadata]][idx_cell_id]
  cell_col_data$cell_type = cell_types
  if (!is.null(plus_to_add_to_metadata)) {
    for (col in plus_to_add_to_metadata) {
      data = metadata[[col]][idx_cell_id]
      cell_col_data[, col] = data
    }
  }
  
  return(cell_col_data)
}

filter_proj_and_add_metadata <- function(proj, nfrags_filter, tss_filter, 
                                         tss_percentile, nfrags_percentile, 
                                         filter_per_cell_type,
                                         dataset, tissue, cell_types, 
                                         min_cells_per_cell_type, 
                                         metadata, 
                                         plus_to_add_to_metadata=NULL,
                                         plus_filters = NULL) {
  cell_col_data = getCellColData(proj)
  
  if (tissue == "all") {
    tissue = "*"
  }
  if (dataset == "all") {
    dataset = "*"
  }
  if (cell_types == "all") {
    cell_types = "*"
  }

  # Filter by Dataset
  print("Filtering by Dataset...")
  dataset_filter = grepl(dataset, cell_col_data[["dataset_per_cell"]])
  proj = proj[dataset_filter]
  cell_col_data = getCellColData(proj)
  print("Done!")
  ##################
  
  # Filter by Tissue
  print("Filtering by tissue...")
  sample_names = unlist(lapply(strsplit(rownames(cell_col_data), split="#"), 
                               "[", 1))
  tissue_filter = grepl(tissue, sample_names)
  proj = proj[tissue_filter]
  cell_col_data = getCellColData(proj)
  print("Done!")
  #################
  
  # Add cell types to cell_col_data
  print("Adding cell types to cell_col_data...")
  cell_col_data = add_cell_types_plus_to_cell_col_data(cell_col_data, metadata, 
                                                  cell_type_col_in_metadata,
                                                  dataset, plus_to_add_to_metadata)
  proj = proj[!is.na(cell_col_data[["cell_type"]])]
  cell_col_data_with_celltypes = cell_col_data[!is.na(cell_col_data[["cell_type"]]), ]
  proj@cellColData = cell_col_data_with_celltypes
  cell_col_data = cell_col_data_with_celltypes
  print("Done!")
  ################################
  
  if (dataset == "Greenleaf_pbmc_bm") {
    pbmc = grepl("PBMC", rownames(cell_col_data))
    cell_col_data[pbmc, "cell_type"] = paste("PBMC", cell_col_data[pbmc, 
                                                                   "cell_type"])
    not_pbmc = !grepl("PBMC", rownames(cell_col_data))
    cell_col_data[not_pbmc, "cell_type"] = paste("Bonemarrow", cell_col_data[not_pbmc, 
                                                                             "cell_type"])
    proj@cellColData = cell_col_data
  } else if (dataset == "Greenleaf_colon") {
    cell_col_data["CellType"] = paste(cell_col_data[["GrossPathology"]], 
                                       cell_col_data[["CellType"]])
    proj@cellColData = cell_col_data
  }
  
  # Filter by cell type
  print("Filtering by cell type...")
  cell_type_filter = grepl(cell_types, cell_col_data[["cell_type"]])
  proj = proj[cell_type_filter]
  cell_col_data = getCellColData(proj)
  print("Done!")
  #####################
  
  # Filter by plus filters
  print("Filtering by additional custom filters...")
  idx = 1
  for (filter in plus_filters) {
    if (!(filter == "NULL")) {
      filter = grepl(filter, cell_col_data[[plus_to_add_to_metadata[idx]]])
      proj = proj[filter]
      cell_col_data = getCellColData(proj)
    }
    idx = idx + 1
  }
  
  #####################
  # Filter by cell number, and TSS/nFrags 
  print("Filtering by TSS/nFrags...")
  counts_per_cell_type = table(cell_col_data[["cell_type"]])
  counts_per_cell_type_filter = counts_per_cell_type >= min_cells_per_cell_type
  cell_types_to_keep = names(counts_per_cell_type)[counts_per_cell_type_filter]
  proj = proj[cell_col_data[["cell_type"]] %in% cell_types_to_keep]
  cell_col_data = getCellColData(proj)
  counts_per_cell_type = table(cell_col_data[["cell_type"]])
  cell_col_data = as.data.frame(cell_col_data)
  if (filter_per_cell_type) {
    cell_col_data = group_by(cell_col_data, cell_type)
  }
  if (!is.null(nfrags_percentile)) {
    temp1 = cell_col_data %>% 
              mutate(throw_away = nFrags < quantile(nFrags, nfrags_percentile))
  } else {
    temp1 = cell_col_data %>% 
      mutate(throw_away = nFrags < nfrags_filter)
  }
  
  temp_frag_filter = temp1[c("throw_away", "cell_type")]
  temp1 = select(temp1, -throw_away)
  temp1 = temp1[!temp_frag_filter[["throw_away"]], ]
  
  if (!is.null(tss_percentile)) {
    temp2 = cell_col_data %>% 
              mutate(throw_away = 
                     TSSEnrichment < quantile(TSSEnrichment, tss_percentile))
  } else {
    temp2 = cell_col_data %>% 
              mutate(throw_away = TSSEnrichment < tss_filter)
  }
  
  temp_tss_filter = temp2[c("throw_away", "cell_type")]
  temp2 = select(temp2, -throw_away)
  temp2 = temp2[!temp_tss_filter[["throw_away"]], ]
  
  temp3 = merge(temp1, temp2)
  
  counts_per_cell_type_after_filtering = table(temp3[["cell_type"]])
  enough_cells = counts_per_cell_type_after_filtering >= min_cells_per_cell_type
  cells_to_filter = names(enough_cells)[enough_cells]
  proj_filter = !(temp_frag_filter["throw_away"] | temp_tss_filter["throw_away"] & 
                 (cell_col_data[["cell_type"]] %in% cells_to_filter))
  
  proj = proj[proj_filter]
  return(proj)
}

reduce_dims <- function(proj, harmonize=F, add_tsne=F, force=F) {
  print("Running iterative LSI")
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(0.2), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=force
  )
  
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  
  proj <- addUMAP(ArchRProj = proj, 
                  reducedDims = "IterativeLSI", 
                  name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine",
                  force=force)
  
  if (add_tsne) {
    proj <- addTSNE(ArchRProj = proj, 
                    reducedDims = "IterativeLSI", 
                    name = "TSNE", nNeighbors = 15,
                    dimsToUse = 2:40,
                    metric = "cosine",
                    force=force)
  }
  
  if (harmonize) {
    proj <- addHarmony(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      name = "Harmony",
      groupBy = "Sample",
      force=force
    )
    
    proj <- addUMAP(ArchRProj = proj, 
                    reducedDims = "Harmony", 
                    name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, 
                    metric = "cosine",
                    force=force)
  }
  proj <- saveArchRProject(ArchRProj = proj,
                           outputDirectory = proj_dir,
                           load = TRUE)
  return(proj)
}

print("Collecting cmd line args")
args = parse_args(OptionParser(option_list=option_list))
# plot_correlation_with_cancer = args$plot_correlation_with_cancer
cores = args$cores
dataset = args$dataset
cluster = args$cluster
cluster_res = args$cluster_res
marker_genes = args$marker_genes
metadata_for_celltype_fn = args$metadata_for_celltype_fn
tissue = args$tissue
nfrags_filter=args$nfrags_filter
tss_filter=args$tss_filter
tss_percentile=args$tss_percentile
nfrags_percentile=args$nfrags_percentile
filter_per_cell_type = args$filter_per_cell_type
cell_types=args$cell_types
plot_custom_column = args$plot_custom_column
sep_for_metadata = args$sep_for_metadata
# cell_name_col_in_metadata = args$cell_name_col_in_metadata
cell_type_col_in_metadata = args$cell_type_col_in_metadata
min_cells_per_cell_type = args$min_cells_per_cell_type
filter_doublets = args$filter_doublets
plot_doublet_scores = args$plot_doublet_scores
save_clusters = args$save_clusters
de_novo_marker_discovery = args$de_novo_marker_discovery
reannotate = args$reannotate
color_embedding_by = args$color_embedding_by
plus_to_add_to_metadata = args$plus_to_add_to_metadata
plus_filters = args$plus_filters
cancer_for_correlation_plot = args$cancer_for_correlation_plot
harmonize = args$harmonize
filter_doublets = args$filter_doublets
doublet_filter = args$doublet_filter
get_metacells = args$get_metacells 

if (!is.null(args$plus_to_add_to_metadata)) {
  plus_to_add_to_metadata = unlist(strsplit(args$plus_to_add_to_metadata, 
                                            split=","))
}

if (!is.null(args$plus_filters)) {
  plus_filters = unlist(strsplit(args$plus_filters, split=","))
}

print("Done collecting cmd line args")

addArchRThreads(threads = cores)
addArchRGenome("hg19")
set.seed(42)
metadata_root = "../../data/metadata"
metadata_filepath = paste(metadata_root, metadata_for_celltype_fn, sep="/")
metadata = read.csv(metadata_filepath, sep=sep_for_metadata)

if (dataset == "Greenleaf_brain") {
  colnames(metadata)[grepl("Iterative.LSI.Clusters", 
                           colnames(metadata))] = "cell_type"
  
  # mapping = list(
  #   c("brain (c0$|c1$|c2$|c5|c6|c7|c13|c14)", "brain GluN"),
  #   c("brain (c3|c4|c12|c16|c20)", "brain IN"), 
  #   c("brain c8", "brain nIPC"),
  #   c("brain c11", "brain early RG"),
  #   c("brain c9", "brain early RG"),
  #   # c("brain (c9|c11)", "brain RG"),
  #   c("brain c10", "brain mGPC"),
  #   c("brain c15", "brain OPC/Oligo"),
  #   c("brain c17", "brain Peric"),
  #   c("brain c19", "brain MG"),
  #   c("brain c21", "brain EC"))
  
  metadata["cell_type"] = gsub("(c0$|c1$|c2$|c5|c6|c7|c13|c14)", 
                               "GluN", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("(c3|c4|c12|c16|c20)", 
                               "IN", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c8", "nIPC", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c9", "late RG", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c11", "early RG", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c10", "mGPC", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c15", "OPC/Oligo", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c17", "Peric", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c19", "MG", metadata[["cell_type"]])
  metadata["cell_type"] = gsub("c21", "EC", metadata[["cell_type"]])
  metadata = metadata[metadata[, "cell_type"] != "c18",]
  # metadata["cell_type"] = gsub("c18", "Unknown", metadata[["cell_type"]])
}

if (dataset == "Shendure" && tissue == "cerebrum") {
  metadata["sample_name"] = gsub("brain", "cerebrum", metadata[["sample_name"]])
}

if (dataset == "Bingren" && tissue == "frontal_cortex") {
  metadata["tissue"] = gsub("Human_brain", "snATAC_frontal_cortex", metadata[["tissue"]])
  metadata["cellID"] = gsub("Human_brain", "snATAC_frontal_cortex", metadata[["cellID"]])
}

setting = paste0("ArchR", "_", "dataset", "_", dataset, "_", "tissue", "_",
                 tissue, "_", "cell_types", "_", cell_types, "_", 
                 "nfrags_filter", "_", nfrags_filter, "_", 
                 "tss_filter", "_", tss_filter, "_", "min_cells_per_cell_type", 
                 "_", min_cells_per_cell_type, "_", "metadata_file", "_", 
                 metadata_for_celltype_fn)

if (!is.null(tss_percentile)) {
  setting = paste0(setting, "_", "tss_percentile", "_", tss_percentile)
}
if (!is.null(nfrags_percentile)) {
  setting = paste0(setting, "_", "nfrags_percentile", "_", nfrags_percentile)
}

if (filter_per_cell_type) {
  setting = paste0(setting, "_", "filter_per_cell_type")
}


idx = 1
for (filter in plus_filters) {
  if (!(filter == "NULL")) {
    setting = paste0(setting, "_", "filter_by_", plus_to_add_to_metadata[idx],
                     "_", filter)
  }
  idx = idx + 1
}

proj_dir = paste("ArchR_projects", setting, sep="/")
setting = gsub("\\|", "_or_", setting)
proj_dir = gsub("\\|", "_or_", proj_dir)

if (dir.exists(proj_dir)) {
  print("Loading existing project")
  proj <- loadArchRProject(proj_dir)
  if (filter_doublets) {
    setting = paste0(setting, "_", "filter_doublets", "_", "doublets_filter", 
                     "_", doublet_filter)
    if (dir.exists(proj_dir)) {
      proj_dir = paste("ArchR_projects", setting, sep="/")
      proj <- loadArchRProject(proj_dir)
    }
    else {
      if (dataset == "Tsankov") {
        cell_col_data = getCellColData(proj)
        idx = cell_col_data[["Clusters_res_0.6"]] == "C4"
        proj = proj[!idx]
        
        # idx_2 = cell_col_data[["DoubletEnrichment"]] >= 10
        # proj = proj[!(idx & idx_2)]
      }
      proj_dir = paste("ArchR_projects", setting, sep="/")
      proj <- saveArchRProject(ArchRProj = proj, 
                               outputDirectory = proj_dir,
                               load = TRUE)
      proj = reduce_dims(proj, harmonize = harmonize, force=T)
    }
  }
} else {
  dir = "ArchR_projects/ArchR_proj"
  print("Loading full ArchR data object")
  ArchR_proj <- loadArchRProject(dir)
  print("Creating new project")
  proj <- filter_proj_and_add_metadata(proj=ArchR_proj, nfrags_filter = nfrags_filter, 
                      tss_filter = tss_filter, tss_percentile = tss_percentile,
                      nfrags_percentile = nfrags_percentile, 
                      filter_per_cell_type = filter_per_cell_type,
                      dataset = dataset, 
                      tissue = tissue, cell_types = cell_types,
                      min_cells_per_cell_type = min_cells_per_cell_type, 
                      metadata = metadata, 
                      plus_to_add_to_metadata = plus_to_add_to_metadata,
                      plus_filters = plus_filters)
  
  print("Saving new project")
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  print("Done saving new project")
  proj = reduce_dims(proj, harmonize, force=T)
}

if (get_metacells) {
  ccd = getCellColData(proj)
  
  ccd["dummy"] = 1
  metaGroupName="dummy"
  # metaGroupName='CellType'
  metaGroup = as.character(ccd[, metaGroupName])
  barcodes = as.character(proj$cellNames)
  metaGroup_df = data.frame(barcode = barcodes, metaGroup = metaGroup)
  k = 500

  KNN = lapply(unique(as.character(ccd[, metaGroupName])), function(x) {
    currentBarcodes = metaGroup_df$barcode[metaGroup_df$metaGroup == x]
    return(as.list(knnGen2(
      ArchRProj = proj,
      reducedDims = 'IterativeLSI',
      overlapCutoff = 0.8,
      cellsToUse = currentBarcodes,
      knnIteration = 10000,
      k = k,
      min_knn = 2,
      seed=10
    )))
  })
  
  # 
  # n = 10
  # knnIteration = currentKnnIteration
  # currentBarcodeslength = length(currentBarcodes)
  # print(paste("Cell type:", x))
  # print(paste("Number of cells:", currentBarcodeslength))
  # currentKnnIteration = floor(currentBarcodeslength / k)
  # currentKnnIteration = 10
  # if (currentKnnIteration == 1) {
  #   currentKnnIteration = 2
  # }
  # k = floor(currentBarcodeslength / n)
  
  # names(KNN) = unique(as.character(ccd[, metaGroupName]))
                      
  # fn = paste(cancer_for_correlation_plot,
  #            "nfrags_filter", nfrags_filter,
  #            "variable_k_n", n, "knnIteration", n*10,
  #            "metacells.Rdata",
  #            sep="_")
  
  fn = paste(dataset,
             "cell_type_independent_nfrags_filter", nfrags_filter,
             "k", k, "knnIteration", "10000",
             "metacells.Rdata",
             sep="_")
  fp = paste("../../data/processed_data", fn, sep="/")
  save(KNN, file=fp)
  
  fn = paste(dataset, "nfrags_filter", nfrags_filter, "embedding.csv", sep="_")
  fp = paste("../../data/processed_data", fn, sep="/")
  
  embedding = "UMAP"
  if (harmonize) {
    embedding = paste0(embedding, "Harmony")
  }
  
  embedding = getEmbedding(proj, embedding)
  write.csv(embedding, fp)
}

# if (plot_correlation_with_cancer) {
#   id = rownames(ccd)
#   ccd = as_tibble(ccd)
#   ccd["id"] = id
#   cancer_correlations_fn = paste(cancer_for_correlation_plot,
#                                  "per_cell_correlations.csv",
#                                  sep="_")
#   cancer_metacorrelations_fn = paste(setting,
#                                      cancer_for_correlation_plot,
#                                      "cell_metacorrelations.csv",
#                                       sep="_")
#   correlations = read.csv(paste("../../data/processed_data", 
#                                 cancer_correlations_fn, sep="/"))
#   metacorrelations = read.csv(paste("../../data/processed_data", 
#                                     cancer_metacorrelations_fn, sep="/"),
#                               row.names = 1)
#   colnames(correlations)[1] = "id"
#   colnames(metacorrelations)[1] = "id"
#   ccd = left_join(ccd, correlations)
#   ccd = left_join(ccd, metacorrelations)
#   ccd = DataFrame(ccd)
#   rownames(ccd) = ccd[, "id"]
#   ccd = ccd[, colnames(ccd) != "id"]
#   proj@cellColData = ccd
#   proj <- saveArchRProject(ArchRProj = proj, 
#                            outputDirectory = proj_dir,
#                            load = TRUE)
#   p <- plotEmbedding(
#     ArchRProj = proj, 
#     colorBy = "cellColData",
#     name = "correlation",
#     embedding = "UMAP",
#     quantCut = c(0.01, 0.95))
#   
#   fn = paste("correlation_with", cancer_for_correlation_plot,
#              setting, sep="_")
#   fn = paste0(fn, ".pdf")
#   print(paste("saving", fn))
#   plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
#   
#   # ecdf_func = ecdf(ccd[["cell_metacorrelation"]])
#   # min_quant = ecdf_func(max(non_zero))
#   
#   ccd[["cell_metacorrelation"]][is.na(ccd[["cell_metacorrelation"]])] = 0
#   non_zero = ccd[["cell_metacorrelation"]][ccd[["cell_metacorrelation"]] < 0]
#   min_val = min(non_zero)
#   max_val = max(non_zero)
#   
#   normalize <- function(x, min_val, max_val) {
#     1 - (x - min_val) / (max_val - min_val)
#   }
#   
#   normalized_values <- normalize(non_zero, min_val, max_val)
#   gradient_func <- colorRamp(c("blue", "red"))
#   gradient_matrix <- gradient_func(normalized_values)
#   gradient_colors <- apply(gradient_matrix, 1, function(row) 
#     rgb(row[1], row[2], row[3], max = 255))
#   
#   ccd[["custom_color"]][ccd[["cell_metacorrelation"]] < 0] = gradient_colors
#   ccd[["custom_color"]][ccd[["custom_color"]] == 0] = "grey"
#   proj@cellColData = ccd
#   
#   p <- plotEmbedding(
#     ArchRProj = proj, 
#     colorBy = "cellColData",
#     name = "custom_color",
#     embedding = "UMAP") + 
#     scale_color_manual(values = unique(ccd[["custom_color"]]))
# 
#   fn = paste("metacorrelation_with", cancer_for_correlation_plot,
#              setting, sep="_")
#   fn = paste0(fn, ".pdf")
#   print(paste("saving", fn))
#   plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
#   
#   p <- plotEmbedding(
#     ArchRProj = proj, 
#     colorBy = "cellColData",
#     name = "nFrags",
#     embedding = "UMAP",
#     quantCut = c(0.01, 0.95))
#   
#   fn = paste("nfrags", cancer_for_correlation_plot,
#              setting, sep="_")
#   fn = paste0(fn, ".pdf")
#   print(paste("saving", fn))
#   plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
# }


if (plot_doublet_scores) {
  proj <- addDoubletScores(
    input = proj,
    useMatrix = "TileMatrix"
  )
  
  # if (dataset == "Tsankov") {
  #   
  # }
  proj = saveArchRProject(ArchRProj = proj, load=TRUE)
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData",
    name = "DoubletEnrichment",
    embedding = "UMAP",
    quantCut = c(0.01, 0.95))
  
  fn = paste("doublets_plot", setting, sep="_")
  fn = paste0(fn, ".pdf")
  print(paste("saving", fn))
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}

if (filter_doublets) {
  ccd = getCellColData(proj)
  if (dataset == "Tsankov") {
    proj = proj[!(ccd[["DoubletEnrichment"]] >= doublet_filter &
                  ccd[["cell_type"]] == "AT2")]
  } else {
    proj = proj[ccd[["DoubletEnrichment"]] < doublet_filter]
  }
  setting = paste(setting, "filter_doublets", "doublet_score", doublet_filter, 
                  sep = "_")
  proj = saveArchRProject(ArchRProj = proj,
                          outputDirectory = setting, 
                          load=TRUE)
  proj = reduce_dims(proj, harmonize, force=T)
}

if (cluster) {
  print(paste0("clustering with clustering resolution = ", cluster_res))
  proj <- addClusters(input = proj,
                      reducedDims = "IterativeLSI",
                      method = "Seurat",
                      name = paste("Clusters", "res", cluster_res, sep="_"),
                      resolution = cluster_res, 
                      force=T)
  if (save_clusters) {
    proj = saveArchRProject(ArchRProj = proj, load=T)
  }
  cell_col_data = getCellColData(proj)
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = paste("Clusters_res", cluster_res, sep="_"),
    embedding = "UMAPHarmony",
    quantCut = c(0.01, 0.95))

  fn = paste("clusters_UMAPs_cluster_res", cluster_res, setting, sep="_")
  # if (filter_doublets) {
  #   fn = paste(fn, "filter_doublets", sep="_")
  # }
  fn = paste0(fn, ".pdf")
  print(paste("saving", fn))
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}

if (reannotate) {
  if (dataset == "Tsankov" && cluster_res==0.6 && cell_types=="Basal") {
    cell_col_data = getCellColData(proj)
    cell_col_data["new_annotation"] = cell_col_data["cell_type"]
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C5", 
                  "new_annotation"] = "Basal.Secretory"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C6" |
                  cell_col_data[["Clusters_res_0.6"]] == "C4", 
                  "new_annotation"] = "Basal.HES2"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C1" |
                  cell_col_data[["Clusters_res_0.6"]] == "C2" |
                  cell_col_data[["Clusters_res_0.6"]] == "C3", 
                  "new_annotation"] = "Basal.TP63+"
    proj = saveArchRProject(ArchRProj = proj, 
                            outputDirectory = setting,
                            load=T)
  }
  if (dataset == "Yang_kidney" && cluster_res==2.5) {
    cell_col_data = getCellColData(proj)
    proj = proj[cell_col_data[["Clusters_res_2.5"]] != "C9" &
      cell_col_data[["Clusters_res_2.5"]] != "C23"]
    setting = paste(setting, "reannotated", sep="_")
    proj = saveArchRProject(ArchRProj = proj, 
                            outputDirectory = proj_dir,
                            load=T)
    fn = paste(dataset, "nfrags_filter", nfrags_filter, "embedding.csv", sep="_")
    fp = paste("../../data/processed_data", fn, sep="/")
    
    embedding = "UMAP"
    if (harmonize) {
      embedding = paste0(embedding, "Harmony")
    }
    
    embedding = getEmbedding(proj, embedding)
    write.csv(embedding, fp)
  }
  if (dataset == "Tsankov" && cluster_res==0.6 && cell_types=="all") {
   cell_col_data = getCellColData(proj)
   # if (filter_doublets) {
   #     idx = cell_col_data[["cell_type"]] == "AT2"
   #     idx_2 = cell_col_data[["DoubletEnrichment"]] >= 10
   #     proj = proj[!(idx & idx_2)]
   # }

   cell_col_data["new_annotation"] = cell_col_data["cell_type"]
   cell_col_data[grepl("Sec-Ciliated", cell_col_data[["new_annotation"]]), 
                  "new_annotation"] = "proximal Ciliated"
   distal_idx = grepl("RPL", rownames(cell_col_data))
   proximal_idx = !distal_idx
   ciliated_idx = grepl("Ciliated", cell_col_data[["cell_type"]])
   secretory_idx = grepl("Secretory", cell_col_data[["cell_type"]])
   cell_col_data[distal_idx & ciliated_idx, "new_annotation"] = "distal_lung Ciliated"
   cell_col_data[proximal_idx & ciliated_idx, "new_annotation"] = "proximal Ciliated"
   cell_col_data[distal_idx & secretory_idx, "new_annotation"] = "distal_lung Secretory"
   cell_col_data[proximal_idx & secretory_idx, "new_annotation"] = "proximal Secretory"
   
    # all_except_idx = !grepl("(proximal|distal) (Ciliated|Secretory)", 
    #                         cell_col_data[["new_annotation"]])
    # cell_col_data[all_except_idx, "new_annotation"] = gsub("(proximal|distal) ", 
    #                                                        "", 
    #                                                        cell_col_data[["new_annotation"]][all_except_idx])
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C3",
                 "new_annotation"] = "Myeloid"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C2", 
                 "new_annotation"] = "Myeloid"
   cell_col_data[cell_col_data[["cell_type"]] == "Fibroblasts", 
                  "new_annotation"] = "Fibroblast.WT1-"
   cell_col_data[cell_col_data[["cell_type"]] == "B_cells", 
                  "new_annotation"] = "T.cells"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C5",
                  "new_annotation"] = "T.cells"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C4", 
                  "new_annotation"] = "B.cells"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C16",
                  "new_annotation"] = "Fibroblast.WT1+"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C17", 
                  "new_annotation"] = "Fibroblast.WT1-"
   cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C1",
                  "new_annotation"] = "Neuronal"
   proj@cellColData = cell_col_data
    
   keep_cells = !(cell_col_data[["new_annotation"]] == "Immune" |
                     cell_col_data[["new_annotation"]] == "Stromal")
   proj = proj[keep_cells, ]
   setting = paste(setting, "reannotated", sep="_")
   proj = saveArchRProject(ArchRProj = proj, 
                            outputDirectory = setting,
                            load=T)
  }
}

if (plot_custom_column) {
  embedding = "UMAP"
  
  if (harmonize) {
    embedding = paste0(embedding, "Harmony")
  }
  
  if (reannotate) {
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "cellColData", 
      name = "new_annotation", 
      embedding = embedding,
      quantCut = c(0.01, 0.95),
      labelMeans = F)
  } else {
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "cellColData", 
      name = color_embedding_by, 
      embedding = embedding,
      quantCut = c(0, 1))
    #quantCut = c(0.01, 0.95))
  }
  if (dataset == "Tsankov" && reannotate) {
    cols <- c("#000075", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
              "#FF0000", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
              "#469990", "#dcbeff", "#9A6324", "#7F00FF", "#800000",
              "#aaffc3", "#808000", "#ffd8b1", "#000000")
    p <- p + 
      scale_color_manual(values = cols,
                         guide = guide_legend(override.aes = 
                                                list(shape = 15)))
  }
  # else {
  # p <- p + 
  #   scale_color_manual(guide = guide_legend(override.aes = 
  #                                             list(shape = 15)))
  # }
  # fn = paste(color_embedding_by, "TSNE", setting, sep="_")
  fn = paste(color_embedding_by, embedding, setting, sep="_")
  
  # if (filter_doublets) {
  #   fn = paste(fn, "filter_doublets", sep="_")
  # }
  if (reannotate) {
    fn = paste(fn, "reannotated", sep="_")
  }
  
  fn = paste0(fn, ".pdf")
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = color_embedding_by, 
    embedding = embedding,
    quantCut = c(0, 1),
    labelMeans=F)
  
  p <- p + 
    ggtitle("") +
          theme(legend.position="none",
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.background = element_rect(fill = '#e0e0e0'))
  
  plotPDF(p, name="fig1.pdf", ArchRProj = proj, addDOC = FALSE)
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = color_embedding_by, 
    embedding = embedding,
    quantCut = c(0, 1),
    labelMeans=F)
  
  p <- p + 
    ggtitle("") +
    theme_classic() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_line(linewidth = 0.1),
          axis.line.y = element_line(linewidth = 0.1))

  
  plotPDF(p, name="fig4A.pdf", ArchRProj = proj, addDOC = FALSE)
  
  if (dataset == "Tsankov") {
    refined_annotation = cell_col_data["new_annotation"]
    write.csv(refined_annotation,
              "../../data/metadata/tsankov_refined_annotation.csv")
  }
}

if (!(is.null(marker_genes))) {
  # impute_weights_dir = paste(proj_dir, "ImputeWeights", sep="/")
  # if (!dir.exists(impute_weights_dir)) {
  proj <- addImputeWeights(proj)
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  # }

  marker_genes = unlist(strsplit(marker_genes, split=","))
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = marker_genes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj),
    quantCut = c(0.01, 0.95),
    height=10,
    width=10,
    baseSize=100
  )
  if (length(marker_genes) == 1) {
    p = list(p)
  }
  p2 <- mapply(function(x, i) {
    x +
      theme_ArchR(baseSize = 5) +
      theme(
        panel.border=element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
        ) +
      guides(fill = guide_colourbar(barwidth = 0.5, barheight = 2))+
      ggtitle(marker_genes[i])
    }, 
    p, seq_along(p), SIMPLIFY = F)
        # plot.title = element_text(hjust = 0.5))
      #   legend.position = "right",
      #   legend.text = element_text(size = 1),
      #   legend.key.size = unit(0.1, "cm")
      # ) +
      # ggtitle(marker_genes[i]) +
      # guides(fill = guide_colourbar(barwidth = 0.5, barheight = 2))
      #guides(color = FALSE, fill = FALSE) 
  # }, p, seq_along(p), SIMPLIFY = F)
  setting = paste0(setting, "_", "marker_genes", "_", paste(marker_genes,
                                                            collapse = "_"))
  path = "../../figures"
  fn = paste0("gene_marker_UMAPs_", setting, ".pdf")
  fp = paste(path, fn, sep="/")
  if (length(marker_genes) > 3) {
    print("Filepath too long")
    print(paste("Saving", fn, "to temp.pdf"))
    fp = paste(path, "temp.pdf", sep="/")
  }
  pdf(fp, width = 10, height = 10)
  do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
  dev.off()
}

if (de_novo_marker_discovery) {
    tryCatch({
      print("Trying to load GeneScoreMatrix")
      gsSE = getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix')
      print("Done loading GeneScoreMatrix")
    }, error = function(e) {
      print("GeneScoreMatrix not found")
      print("Creating GeneScoreMatrix")
      proj = addGeneScoreMatrix(proj)
      proj = saveArchRProject(ArchRProj = proj, 
                               outputDirectory = proj_dir,
                               load = TRUE)
    })
    gsSE = gsSE[, proj$cellNames]
    metaGroupName = paste("Clusters", "res", cluster_res, sep="_")
    if (!file.exists(paste0(proj_dir, 'DAG_', metaGroupName, '.rds'))) {
      DAG_list = getMarkerFeatures(ArchRProj = proj,
                                   testMethod = "wilcoxon",
                                   binarize = FALSE,
                                   useMatrix = "GeneScoreMatrix",
                                   groupBy = metaGroupName)
      listnames = colnames(DAG_list)
      DAG_list = lapply(1:ncol(DAG_list), function(x) 
        {
          df = DAG_list[, x]
          df = do.call(cbind, (assays(df)))
          colnames(df) = names(assays(DAG_list))
          df$gene = rowData(DAG_list)$name
          return(df)
        })
      names(DAG_list) = listnames
      saveRDS(DAG_list, paste0(proj_dir, '/DAG_', metaGroupName,'.rds'))
  } else {
    DAG_list = readRDS(paste0(proj_dir,'/DAG_', metaGroupName,'.rds'))
  }
  FDR_threshold = 1e-2
  lfc_threshold = 2
  top_genes = 50
  
  DAG_top_list = DAG_list[sapply(DAG_list, function(x) {
                                              nrow(x[x$FDR < FDR_threshold & 
                                                       abs(x$Log2FC) > 
                                                       lfc_threshold,]) > 0
                                           })]
                                  
  DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
    res = DAG_top_list[[x]]
    res = na.omit(res)
    res = res[res$FDR < FDR_threshold, ]
    res = res[order(res$FDR), ]
    res = res[abs(res$Log2FC) > lfc_threshold,]
    res$comparison = names(DAG_top_list)[x]
    if (nrow(res) < top_genes) { 
      res 
      } else { 
        head(res, top_genes) 
      }
  })
  DAG_df = Reduce(rbind, DAG_top_list)
  
  gsMat = assays(gsSE)[[1]]
  rownames(gsMat) = rowData(gsSE)$name
  gsMat_mg = gsMat[rownames(gsMat) %in% DAG_df$gene, ]
  gsMat_mg = as.data.frame(t(gsMat_mg))
  gsMat_mg$metaGroup = as.character(proj@cellColData[,metaGroupName])
  gsMat_mg = aggregate(. ~ metaGroup, gsMat_mg, mean)
  rownames(gsMat_mg) = gsMat_mg[, 1]
  gsMat_mg = gsMat_mg[, -1]
  gsMat_mg = gsMat_mg[names(table(proj@cellColData[, metaGroupName])[table 
                                    (proj@cellColData[, metaGroupName]) > 50]),]
  DAG_hm = Heatmap(t(scale(gsMat_mg)), 
                    row_labels = colnames(gsMat_mg),
                    column_title = paste('top', top_genes),
                    clustering_distance_columns = 'euclidean',
                    clustering_distance_rows = 'euclidean',
                    cluster_rows = T,
                    #col = pals_heatmap[[5]],
                    cluster_columns=T,#col = pals_heatmap[[1]],
                    row_names_gp = gpar(fontsize = 2),
                    column_names_gp = gpar(fontsize = 4),
                    rect_gp = gpar(col = "white", lwd = .5),
                    border=TRUE
                    #right_annotation = motif_ha
  )
  
  DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', 
                                column_title_gp = gpar(fontsize = 13)))
  pdf(paste0(proj_dir, '/Plots/DAG_clusters_', metaGroupName, '_heatmaps.pdf'), 
       width = 4, height = 10)
  DAG_hm
  dev.off()
}

########## Post visual inspection ##########
############################################
############################################
map <- function(x, mapping) mapping[[x]]

if (dataset == "Tsankov" && tissue == "all" && nfrags_filter == 1000 && 
  tss_filter == 4 && cell_types == "Basal") {
  mapping = c("C1" = "Basal,Sec.SOX4", "C2" = "Basal.HES2", 
              "C3" = "Basal.TP63,SOX2", "C4" = "Basal.TP63",
              "C5" = "Basal.FOXA1")
  
  mapped_cells = unlist(lapply(cell_col_data[["Clusters"]], map, mapping))
  new_annotation_df = data.frame(Sample = rownames(cell_col_data), 
                                 celltypes = mapped_cells)
  
  dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj"
  ArchR_proj <- loadArchRProject(dir)
  cell_col_data = getCellColData(ArchR_proj)
  cell_col_data = add_cell_types_to_cell_col_data(cell_col_data, metadata,
                                                  cell_type_col_in_metadata,
                                                  dataset)
  ArchR_proj@cellColData = cell_col_data
  # FIX 
  # proj <- filter_proj(ArchR_proj, nfrags_filter, tss_filter, tss_percentile,
  #                     dataset, "IC", "all", metadata)
  cell_col_data = getCellColData(proj)
  non_basal_proximal = cell_col_data[cell_col_data[["cell_type"]] != "Basal", ]
  non_basal_proximal = as.data.frame(non_basal_proximal["cell_type"])
  colnames(non_basal_proximal) = "celltypes"
  rownames(new_annotation_df) = new_annotation_df$Sample
  new_annotation_df$Sample= NULL
  
  new_annotation_df = rbind(new_annotation_df, non_basal_proximal)
  annotation_filename = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata", 
                        paste0(setting, "_annotation.csv"), sep="/")
  write.csv(new_annotation_df, annotation_filename)
  
  #-- Call peaks --#
  cellGroup = "Clusters"
  proj = addGroupCoverages(
    ArchRProj = proj,
    groupBy = cellGroup,
    force = TRUE,
    minCells= 20,
    maxCells = 500,
    minReplicates = 2,
    sampleRatio = 0.8,
    useLabels = TRUE)
  
  proj = addReproduciblePeakSet(
    proj,
    groupBy= cellGroup,
    peakMethod = 'Macs2',
    reproducibility = "2",
    maxPeaks = 500000,
    minCells=20, force=TRUE)
  
  proj = addPeakMatrix(proj)
  
  proj <- addMotifAnnotations(ArchRProj = proj, 
                              motifSet = "cisbp", name = "Motif",
                              force=TRUE)
  proj <- addBgdPeaks(proj)
  proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  saveArchRProject(ArchRProj = proj)
  
  markers = c("TP63", "SOX2", "HES2", "FOXA1", "SOX4", "NKX21")
  markers_p = list() 
  for (i in markers) {
    markerMotifs = getFeatures(proj, select = paste0(i, "_"), 
                               useMatrix = "MotifMatrix")
    markerMotifs2 = grep("z:", markerMotifs, value = TRUE)
    markers_p[[i]] = plotEmbedding(ArchRProj = proj, 
                                    colorBy = "MotifMatrix", 
                                    name = markerMotifs2, 
                                    embedding = "UMAP", 
                                    size=1, 
                                    baseSize=0, 
                                    imputeWeights = getImputeWeights(proj), 
                                    plotAs="points")
  }
  markers_p2 = lapply(1:length(markers_p), function(x) { 
    markers_p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) + 
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
      theme(axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank() 
      )})
  png(paste0(proj_dir, "/Plots/variable_motifs.png"), 
      width=800, height=1000)
  wrap_plots(markers_p2, ncol=6)
  dev.off()
}




# out = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis", 
#             "ArchR_per_dataset", dataset, tissue, sep="/")
# dir.create(out, recursive=T)
# saveArchRProject(ArchRProj = proj, outputDirectory=out)

