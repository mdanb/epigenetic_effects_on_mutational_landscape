library(ArchR)
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(dplyr)
library(ComplexHeatmap)

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character", default="all"),
  make_option("--metadata_for_celltype_fn", type="character"),
  make_option("--sep_for_metadata", type="character", default=","),
  make_option("--cell_type_col_in_metadata", type="character"),
  # make_option("--cell_name_col_in_metadata", type="character"),
  make_option("--cluster", action="store_true", default=F),
  make_option("--cluster_res", type="double"),
  make_option("--plot_cell_types", action="store_true", default=F),
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
  make_option("--de_novo_marker_discovery", action="store_true", default=FALSE)
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
#                       "--marker_genes=SCGB3A1,SCGB3A2,TP63,KRT5")
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
#                       "--filter_doublets",
#                       "--marker_genes=CEBPE,GFI1,IRF1,GNL2,ELANE,PPARG,BHLHE41,EGR2,FABP4,HBA1,MARCO,MME,F13A1,SLC40A1,SELP,STAB1,FOLR2,CD1C,CD1A,FCER1A,HLA-DQA1,CLEC10A,PKIB,CLEC9A,GCSAM,BATF3,WDFY4,LILRA4,GZMB,IL3RX,LAMP3,CCR7,FSCN1,CCL22,MARCKSL1,EBI3,IDO1,S100A8,S100A9,S100A12,FCN1,CD14,VCAN,PTX3,LRG1,LYPD2,LST1,LILRB2,FCGR3A,NKG7,NCR1,SPON2,KLRD1,GNLY,KLRC1,KLRF1,FGFBP2,KIT,TPSAB1,CPA3,CTSG,BACH2,ITGA2B,GP1BA,VWF,FLI1,MS4A1,SDC1,KLK1,EBF1,CD79A,CD79B,PAX5,VPREB3,POU2AF1,MZB1,XBP1,CD3D,LEF1,TCF7,TCF3,IL7R,CD8A,CD4,BCL11B,FOXP3,IL2RA,TNFRSF4,TNFRSF18,CTLA4,EPCAM,CLU,CLDN18,HNF1B,PECAM1,KDR,PLVAP,MCAM,COL1A1,COL1A2,SPARCL1,PDGFRA,NEUROG1,OLIG2,RFX2,NOTO,POU2F1,RFX5,C1QA,C1QC,CD68,CD163,THBS1,FN1,C5AR2")
# )
args = parse_args(OptionParser(option_list=option_list), args=
                    c("--cores=8",
                      "--dataset=Tsankov",
                      "--metadata_for_celltype_fn=combined_distal_proximal.csv",
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
                      "--reannotate"
))

add_cell_types_to_cell_col_data <- function(cell_col_data, metadata,
                                            cell_type_col_in_orig_metadata, 
                                            dataset) {
  if (dataset == "Shendure") {
    to_match = paste(metadata[["sample_name"]], metadata[["cell"]], sep="#")
    rownames_archr = rownames(cell_col_data)
  }
  else if (dataset == "Bingren") {
    metadata = metadata[metadata["Life.stage"] == "Adult", ]
    metadata["cellID"] = gsub("\\+", "#", metadata[["cellID"]])
    idx = str_locate(metadata[["cellID"]], "#")[, 1] 
    to_match = unname(mapply(function(s, i, ins) {
      paste0(substring(s, 1, i-1), ins, substring(s, i))
    }, s = metadata[["cellID"]], i = idx-1, ins = "rep"))
    rownames_archr = unlist(lapply(lapply(strsplit(rownames(cell_col_data), split="_"), 
                         "[", 2:4), paste, collapse="_"))
  }
  else if (dataset == "Tsankov") {
    rownames_archr = rownames(cell_col_data)
    to_match = metadata[["X"]]
  }
  idx_cell_id = match(rownames_archr, to_match)
  cell_types = metadata[[cell_type_col_in_orig_metadata]][idx_cell_id]
  cell_col_data$cell_type = cell_types
  return(cell_col_data)
}

filter_proj <- function(proj, nfrags_filter, tss_filter, tss_percentile,
                        nfrags_percentile, filter_per_cell_type,
                        dataset, tissue, cell_types, min_cells_per_cell_type, 
                        metadata) {
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
  dataset_filter = grepl(dataset, cell_col_data[["dataset_per_cell"]])
  proj = proj[dataset_filter]
  cell_col_data = getCellColData(proj)
  ##################
  
  # Filter by Tissue
  sample_names = unlist(lapply(strsplit(rownames(cell_col_data), split="#"), 
                               "[", 1))
  tissue_filter = grepl(tissue, sample_names)
  proj = proj[tissue_filter]
  cell_col_data = getCellColData(proj)
  #################
  
  # Add cell types to cell_col_data
  cell_col_data = add_cell_types_to_cell_col_data(cell_col_data, metadata, 
                                                  cell_type_col_in_metadata,
                                                  dataset)
  proj = proj[!is.na(cell_col_data[["cell_type"]])]
  cell_col_data_with_celltypes = cell_col_data[!is.na(cell_col_data[["cell_type"]]), ]
  proj@cellColData = cell_col_data_with_celltypes
  cell_col_data = cell_col_data_with_celltypes
  ################################
  
  # Filter by cell type
  cell_type_filter = grepl(cell_types, cell_col_data[["cell_type"]])
  proj = proj[cell_type_filter]
  cell_col_data = getCellColData(proj)
  #####################
  
  # Filter by cell number, and TSS/nFrags 
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
              mutate(throw_away = nFrags < tss_filter)
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

reduce_dims <- function(proj, force=F) {
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
  
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  return(proj)
}

print("Collecting cmd line args")
args = parse_args(OptionParser(option_list=option_list))
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
plot_cell_types = args$plot_cell_types
sep_for_metadata = args$sep_for_metadata
# cell_name_col_in_metadata = args$cell_name_col_in_metadata
cell_type_col_in_metadata = args$cell_type_col_in_metadata
min_cells_per_cell_type = args$min_cells_per_cell_type
filter_doublets = args$filter_doublets
plot_doublet_scores = args$plot_doublet_scores
save_clusters = args$save_clusters
de_novo_marker_discovery = args$de_novo_marker_discovery
reannotate = args$reannotate
print("Done collecting cmd line args")

addArchRThreads(threads = cores)
addArchRGenome("hg19")
set.seed(42)
metadata_root = "../../data/metadata"
metadata_filepath = paste(metadata_root, metadata_for_celltype_fn, sep="/")
metadata = read.csv(metadata_filepath, sep=sep_for_metadata)

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

proj_dir = paste("ArchR_projects", setting, sep="/")

if (dir.exists(proj_dir)) {
  print("Loading existing project")
  proj <- loadArchRProject(proj_dir)
  if (filter_doublets) {
    setting = paste0(setting, "_", "filter_doublets")
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
      proj = reduce_dims(proj, force=T)
    }
  }
} else {
  dir = "ArchR_projects/ArchR_proj/"
  print("Loading full ArchR data object")
  ArchR_proj <- loadArchRProject(dir)
  print("Creating new project")
  proj <- filter_proj(proj = ArchR_proj, nfrags_filter, tss_filter, tss_percentile, 
                      nfrags_percentile, filter_per_cell_type, 
                      dataset, tissue, cell_types,
                      min_cells_per_cell_type, metadata)
  
  print("Saving new project")
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  print("Done saving new project")
  proj = reduce_dims(proj)
}

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
    embedding = "UMAP",
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
  if (dataset == "Tsankov" && cluster_res==0.6) {
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
    cell_col_data[cell_col_data[["cell_type"]] == "Fibroblasts", 
                  "new_annotation"] = "Fibroblast.WT1-"
    cell_col_data[cell_col_data[["cell_type"]] == "B_cells", 
                  "new_annotation"] = "B.cells"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C18", 
                  "new_annotation"] = "T.cells"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C17", 
                  "new_annotation"] = "B.cells"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C5", 
                  "new_annotation"] = "Fibroblast.WT1+"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C6", 
                  "new_annotation"] = "Fibroblast.WT1-"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C9", 
                  "new_annotation"] = "Myeloid"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C12", 
                  "new_annotation"] = "Myeloid"
    cell_col_data[cell_col_data[["Clusters_res_0.6"]] == "C1", 
                  "new_annotation"] = "Neuronal"
    proj@cellColData = cell_col_data
    
    keep_cells = !(cell_col_data[["new_annotation"]] == "Immune" |
                     cell_col_data[["new_annotation"]] == "Stromal")
    proj = proj[keep_cells, ]
    proj = saveArchRProject(ArchRProj = proj, load=T)
  }
}

if (plot_cell_types) {
  if (reannotate) {
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "cellColData", 
      name = "new_annotation", 
      embedding = "UMAP",
      quantCut = c(0.01, 0.95))
  } else {
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "cellColData", 
      name = "cell_type", 
      embedding = "UMAP",
      quantCut = c(0.01, 0.95))
  }
  cols <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
            "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
            "#469990", "#dcbeff", "#9A6324", "#7F00FF", "#800000",
            "#aaffc3", "#808000", "#ffd8b1", "#000075", "#000000")
  p <- p + 
    scale_color_manual(values = cols,
                       guide = guide_legend(override.aes = 
                                              list(shape = 15)))  
  fn = paste("cell_type_UMAP", setting, sep="_")
  # if (filter_doublets) {
  #   fn = paste(fn, "filter_doublets", sep="_")
  # }
  if (reannotate) {
    fn = paste(fn, "reannotated", sep="_")
  }
  fn = paste0(fn, ".pdf")
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
  
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
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 10) +
      theme(plot.margin = unit(c(0.5, 0, 0.5, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line=element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank()
      ) 
  })
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
  pdf(fp, width = 20, height = 20)
  do.call(cowplot::plot_grid, c(list(ncol = 5), p2))
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
  markers_p2 = lapply(1:length(markers_p), function(x){ 
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

