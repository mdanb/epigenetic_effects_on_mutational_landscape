library(ArchR)
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
# make_option("--function_number_to_process_celltype", type="integer"),

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character", default="all"),
  make_option("--metadata_for_celltype_fn", type="character"),
  make_option("--sep_for_metadata", type="character", default=","),
  make_option("--cell_type_col_in_metadata", type="character"),
  make_option("--cell_name_col_in_metadata", type="character"),
  make_option("--cluster", action="store_true", default=F),
  make_option("--cluster_res", type="double", default=1.2),
  make_option("--plot_cell_types", action="store_true", default=F),
  make_option("--tissue", type="character", default="all"),
  make_option("--nfrags_filter", type="integer", default=1),
  make_option("--tss_filter", type="integer", default=0),
  make_option("--tss_percentile", type="double", default=NULL),
  make_option("--nfrags_percentile", type="double", default=NULL),
  make_option("--cell_types", type="character", default="all"),
  make_option("--marker_genes", type="character", default=NULL),
  make_option("--min_cells_per_cell_type", type="integer", default=0)
)

# # Tsankov combined distal, proximal with epithelial gene markers
# args = parse_args(OptionParser(option_list=option_list), args=
#                   c("--cores=4",
#                     "--dataset=Tsankov",
#                     "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                     "--sep_for_metadata=,",
#                     "--cell_type_col_in_metadata=celltypes",
#                     "--cell_name_col_in_metadata=X",
#                     "--column_to_color_by=NULL",
#                     "--tissue=all",
#                     "--nfrags_filter=1000",
#                     "--tss_filter=4",
#                     "--cell_types=all",
#                     "--marker_genes=KRT15,KRT17,KRT5,S100A2,EPCAM,KRT4,KRT13,TP63,SOX2,HES2,FOXA1,SOX4,NKX2-1,SCGB1A1,SCGB3A1,SCGB3A2,MUC5B"))
# 
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=Tsankov_fibro-fibro+C12+fibro+C14.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cell_name_col_in_metadata=X",
#                       "--column_to_color_by=NULL",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--marker_genes=WT1"))
# 
# # Basal
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cell_name_col_in_metadata=X",
#                       "--column_to_color_by=NULL",
#                       "--tissue=all",
#                       "--nfrags_filter=1000",
#                       "--tss_filter=4",
#                       "--cell_types=Basal",
#                       "--marker_genes=KRT15,KRT17,KRT5,S100A2,EPCAM,KRT4,KRT13,TP63,SOX2,HES2,FOXA1,SOX4,NKX2-1,SCGB1A1,SCGB3A1,SCGB3A2,MUC5B"))
#                       
# # Shendure cerebellum
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Shendure",
#                       "--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--cell_name_col_in_metadata=cell",
#                       "--column_to_color_by=cell_type",
#                       "--tissue=cerebellum",
#                       "--nfrags_filter=1000",
#                       "--tss_percentile=0.25",
#                       "--tss_filter=4",
#                       "--cell_types=all")
#                   )
# 
# # Shendure cerebrum
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Shendure",
#                       "--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--cell_name_col_in_metadata=cell",
#                       "--column_to_color_by=cell_type",
#                       "--tissue=cerebrum",
#                       "--nfrags_filter=1000",
#                       "--tss_percentile=0.2",
#                       "--nfrags_percentile=0.2",
#                       "--tss_filter=4",
#                       "--cell_types=all")
#                   )
# ,
#                       "--marker_genes=AQP4,OLIG1,OLIG2"
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Shendure",
#                       "--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--cell_name_col_in_metadata=cell",
#                       "--column_to_color_by=cell_type",
#                       "--tissue=cerebrum",
#                       "--nfrags_filter=1",
#                       "--tss_filter=1",
#                       "--cell_types=all")
# )
# 
# # Bing Ren frontal cortex
# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=4",
#                       "--dataset=Bingren",
#                       "--metadata_for_celltype_fn=GSE184462_metadata.tsv",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell.type",
#                       "--cell_name_col_in_metadata=cellID",
#                       "--column_to_color_by=cell.type",
#                       "--tissue=frontal_cortex",
#                       "--nfrags_filter=1000",
#                       "--tss_filter=4",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=50")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Bingren",
#                       "--metadata_for_celltype_fn=bingren_remove_same_celltype_indexing.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=cell.type",
#                       "--cell_name_col_in_metadata=cellID",
#                       "--column_to_color_by=cell.type",
#                       "--tissue=frontal_cortex",
#                       "--cell_types=all",
#                       "--marker_genes=OLIG1,OLIG2,AQP4,RBFOX3,EOMES",
#                       "--min_cells_per_cell_type=100")
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cell_name_col_in_metadata=X",
#                       "--column_to_color_by=NULL",
#                       "--tissue=all",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--cell_types=all",
#                       "--marker_genes=OLIG1,OLIG2,AQP4,RBFOX3,EOMES",
#                       "--min_cells_per_cell_type=100")
# )

# args = parse_args(OptionParser(option_list=option_list))



# func_1 <- function(original_metadata, cell_col_data, cell_type_col_in_metadata,
#                    cell_name_col_in_metadata) {
#   cell_id = strsplit(original_metadata[[cell_name_col_in_metadata]], "#")
#   cell_id = lapply(cell_id, "[", 2)
#   cell_id = substr(cell_id, 1, 16)
#   df[cell_id] = cell_id
#   return(df[, c(cell_name_col_in_metadata, cell_type_col_in_metadata)])
# }

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Bingren",
#                       "--metadata_for_celltype_fn=GSE184462_metadata.tsv",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell.type",
#                       "--cell_name_col_in_metadata=cellID",
#                       "--tissue=stomach",
#                       "--nfrags_percentile=0.2",
#                       "--tss_percentile=0.2",
#                       "--cell_types=all",
#                       "--marker_genes=AGR2,CLCA1,KLF4,MUC2,MUC5B,SPDEF,TFF3,ATP4A,FUT2,MUC6,REG4,TFF1",
#                       "--min_cells_per_cell_type=200"
#                     )
# )

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Shendure",
#                       "--metadata_for_celltype_fn=GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
#                       "--sep_for_metadata=\t",
#                       "--cell_type_col_in_metadata=cell_type",
#                       "--cell_name_col_in_metadata=cell",
#                       "--column_to_color_by=cell_type",
#                       "--tissue=stomach",
#                       "--cell_types=all",
#                       "--marker_genes=TFF1,MUC2,TFF3,ATP4A,MUC5B,CLCA1,KLF4,MUC6,FUT2,REG4,AGR2,SPDEF",
#                       "--min_cells_per_cell_type=100")
# )


add_cell_types_to_cell_col_data <- function(cell_col_data, metadata,
                                            cell_type_col_in_orig_metadata, 
                                            dataset) {
  if (dataset == "Shendure") {
    to_match = paste(metadata[["sample_name"]], metadata[["cell"]], sep="#")
    rownames_archr = rownames(cell_col_data)
    
    # cell_id = original_metadata[[cell_name_col_in_orig_metadata]]
    # archr_name = paste(original_metadata[["sample_name"]],
    #                    original_metadata[[cell_name_col_in_orig_metadata]],
    #                    sep="#")
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
    # cell_id = unlist(lapply(strsplit(metadata[["cellID"]], "\\+"), 
    #                           "[", 2))
    # sample_to_match_archr_in_metadata = unlist(lapply(lapply(strsplit(cell_col_data[["Sample"]], 
    #                                         split = "_"), "[", 2:3), paste,
    #                                     collapse = "_"))
  }
  else if (dataset == "Tsankov") {
    rownames_archr = rownames(cell_col_data)
    to_match = metadata[["X"]]
  }
  # else if (dataset == "Bingren") {
  #   temp = original_metadata[[cell_name_col_in_orig_metadata]]
  #   index <- unlist(gregexpr("\\+", temp))
  #   part1 <- substring(temp, 1, index-4)
  #   part2 <- substring(temp, index-3, index-3)
  #   part3 <- substring(temp, index+1, nchar(temp))
  #   
  #   archr_name = paste0("GSM5589414_UMB4540_", part1, "rep", part2, "#", part3)
  # }
  # else {
  #   archr_name = unlist(original_metadata[[cell_name_col_in_orig_metadata]])
  # }
  # archr_cells = unlist(lapply(strsplit(rownames(cell_col_data), "#"), "[", 2))
  
  idx_cell_id = match(rownames_archr, to_match)
  cell_types = metadata[[cell_type_col_in_orig_metadata]][idx_cell_id]
  # not_na_idx = !is.na(idx)
  # cell_types = cell_types[not_na_idx]
  # cell_col_data[not_na_idx, "cell_type"] = cell_types
  cell_col_data$cell_type = cell_types
  return(cell_col_data)
}

filter_proj <- function(proj, nfrags_filter, tss_filter, tss_percentile, nfrags_percentile,
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

  if (!is.null(tss_percentile)) {
    cutoff <- quantile(cell_col_data[["TSSEnrichment"]], tss_percentile)
    tss_filter = cell_col_data[["TSSEnrichment"]] >= cutoff
  } else {
    tss_filter = cell_col_data[["TSSEnrichment"]] >= tss_filter
  }
  
  if (!is.null(nfrags_percentile)) {
    cutoff <- quantile(cell_col_data[["nFrags"]], tss_percentile)
    frags_filter = cell_col_data[["nFrags"]] >= cutoff
  } else {
    frags_filter = cell_col_data[["nFrags"]] >= nfrags_filter
  }
  
  dataset_filter = grepl(dataset, cell_col_data[["dataset_per_cell"]])
  
  sample_names = unlist(lapply(strsplit(rownames(cell_col_data), split="#"), 
                               "[", 1))
  tissue_filter = grepl(tissue, sample_names)
  proj = proj[frags_filter & tss_filter & dataset_filter & tissue_filter]
  cell_col_data = getCellColData(proj)
  # cell_col_data = cell_col_data[grepl("stomach_SM-CHLWL", rownames(cell_col_data)), ]
  
  cell_col_data = add_cell_types_to_cell_col_data(cell_col_data, metadata, 
                                                  cell_type_col_in_metadata,
                                                  dataset)
  proj = proj[!is.na(cell_col_data[["cell_type"]])]
  proj@cellColData = cell_col_data[!is.na(cell_col_data[["cell_type"]]), ]
  cell_col_data = getCellColData(proj)
  cell_type_filter = grepl(cell_types, cell_col_data[["cell_type"]])
  proj = proj[cell_type_filter]
  counts_per_cell_type = table(cell_col_data[["cell_type"]])
  counts_per_cell_type_filter = counts_per_cell_type >= min_cells_per_cell_type
  cell_types_to_keep = names(counts_per_cell_type)[counts_per_cell_type_filter]
  proj = proj[cell_col_data[["cell_type"]] %in% cell_types_to_keep]
  return(proj)
}

# args = parse_args(OptionParser(option_list=option_list))

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
cell_types=args$cell_types
plot_cell_types = args$plot_cell_types
sep_for_metadata = args$sep_for_metadata
# column_to_color_by = args$column_to_color_by

cell_name_col_in_metadata = args$cell_name_col_in_metadata
cell_type_col_in_metadata = args$cell_type_col_in_metadata
min_cells_per_cell_type = args$min_cells_per_cell_type
print("Done collecting cmd line args")

addArchRThreads(threads = cores)
addArchRGenome("hg19")
set.seed(42)
# arrow_files_path = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/arrow",
#                           dataset, sep="/")

# arrow_files = list.files(arrow_files_path, full.names=T, pattern="arrow")
# markers = c('WT1','ITLN1','COL1A1','PECAM1','LYZ','CD3D','MSLN','KRT18','KRT5','VIM')
metadata_root = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata"
metadata_filepath = paste(metadata_root, metadata_for_celltype_fn, sep="/")
metadata = read.csv(metadata_filepath, sep=sep_for_metadata)

if (dataset == "Shendure" && tissue == "cerebrum") {
  metadata["sample_name"] = gsub("brain", "cerebrum", metadata[["sample_name"]])
}

if (dataset == "Bingren" && tissue == "frontal_cortex") {
  metadata["tissue"] = gsub("Human_brain", "snATAC_frontal_cortex", metadata[["tissue"]])
  metadata["cellID"] = gsub("Human_brain", "snATAC_frontal_cortex", metadata[["cellID"]])
}

root = paste0("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis")
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

proj_dir = paste(root, setting, sep="/")

if (dir.exists(proj_dir)) {
  print("Loading existing project")
  proj <- loadArchRProject(proj_dir)
  print("Done loading existing project")
} else {
  print("Creating new project")
  dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj"
  ArchR_proj <- loadArchRProject(dir)
  proj <- filter_proj(proj = ArchR_proj, nfrags_filter, tss_filter, tss_percentile, 
                      nfrags_percentile, dataset, tissue, cell_types,
                      min_cells_per_cell_type, metadata)
  
  print("Saving new project")
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  print("Done saving new project")
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
    dimsToUse = 1:30
  )
  
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  
  proj <- addUMAP(ArchRProj = proj, 
                  reducedDims = "IterativeLSI", 
                  name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                  metric = "cosine")
  
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
}
# tryCatch({
# getMatrixFromProject(ArchR_proj_subset, useMatrix = "GeneScoreMatrix")
# }, error = function(err) {
# proj <- addGeneScoreMatrix(ArchR_proj_subset, force=T)
# })

# tryCatch({
  # getMatrixFromProject(ArchR_proj_subset, useMatrix = "TileMatrix")
# }, error = function(err) {
# ArchR_proj_subset <- addTileMatrix(ArchR_proj_subset, force=T)
# })

if (!(is.null(marker_genes))) {
  proj <- addImputeWeights(proj)
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  marker_genes = unlist(strsplit(marker_genes, split=","))
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = marker_genes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj),
    quantCut = c(0.01, 0.95),
    height=10,
    width=10
  )
  proj <- saveArchRProject(ArchRProj = proj, 
                           outputDirectory = proj_dir,
                           load = TRUE)
  if (length(marker_genes) == 1) {
    p = list(p)
  }
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 3) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
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
  path = "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/figures"
  fn = paste0("gene_marker_UMAPs_", setting, ".pdf")
  fp = paste(path, fn, sep="/")
  if (length(marker_genes) > 3) {
     fp = paste(path, "temp.pdf", sep="/")
  }
  print("Filepath too long")
  print(paste("Saving", fp, "to temp.pdf"))
  pdf(fp)
  do.call(cowplot::plot_grid, c(list(ncol = 5), p2))
  dev.off()
}

if (cluster) {
# if (!("Clusters" %in% colnames(cell_col_data))) {
  print(paste0("clustering with clustering resolution = ", cluster_res))
  proj <- addClusters(input = proj,
                      reducedDims = "IterativeLSI",
                      method = "Seurat",
                      name = "Clusters",
                      resolution = cluster_res, 
                      force=T)
  cell_col_data = getCellColData(proj)
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "cellColData", 
    name = "Clusters",
    embedding = "UMAP",
    quantCut = c(0.01, 0.95))
  fn = paste0("clusters_UMAPs", "_", "dataset", "_", dataset, "_", "tissue", "_",
              tissue, "_", "cell_types", "_", cell_types, "_", "cluster_res", "_", 
              cluster_res, "_", "nfrags_filter", "_", nfrags_filter, "_",  
              "tss_filter", "_", tss_filter)
  if (!is.null(tss_percentile)) {
    fn = paste0(fn, "_tss_percentile_", tss_percentile)
  } 
  if (!is.null(nfrags_percentile)) {
    fn = paste0(fn, "_nfrags_percentile_", nfrags_percentile)
  } 
  if (!is.null(min_cells_per_cell_type)) {
    fn = paste0(fn, "_min_cells_per_cell_type_", min_cells_per_cell_type)
  } 
  fn = paste0(fn, ".pdf")
  print(paste("saving", fn))
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}
  # if (dataset == "Tsankov") {
  #   write.csv(cell_col_data["Clusters"], 
  #   "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/metadata_Tsankov_refined_fibroblasts.csv")
  # }
if (plot_cell_types) {
  p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "cellColData", 
        name = "cell_type", 
        embedding = "UMAP",
        quantCut = c(0.01, 0.95))
  fn = paste0("cell_type_UMAP", "_", "dataset", "_", dataset, "_", "tissue", "_",
              tissue, "_", "cell_types", "_", cell_types, "_", 
              "nfrags_filter", "_", nfrags_filter, "_",  
              "tss_filter", "_", tss_filter)
  
  if (!is.null(tss_percentile)) {
    fn = paste0(fn, "_tss_percentile_", tss_percentile)
  } 
  if (!is.null(nfrags_percentile)) {
    fn = paste0(fn, "_nfrags_percentile_", nfrags_percentile)
  } 
  if (!is.null(min_cells_per_cell_type)) {
    fn = paste0(fn, "_min_cells_per_cell_type_", min_cells_per_cell_type)
  } 
  fn = paste0(fn, ".pdf")
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}

# Post visual inspection
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
                                                  cell_name_col_in_metadata,
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
    markerMotifs = getFeatures (proj, select = paste0(i, "_"), 
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
  # dev = getMatrixFromProject(proj, useMatrix = "MotifMatrix")
  # TF = 'z:PITX2_504'
  # p <- plotGroups(ArchRProj = proj,
  #                 groupBy = "Clusters", 
  #                 colorBy = "MotifMatrix",
  #                 name = TF,
  #                 seqnames ='z',
  #                 imputeWeights = getImputeWeights(proj)
  #                )
}

saveArchRProject(ArchRProj = proj)



# out = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis", 
#             "ArchR_per_dataset", dataset, tissue, sep="/")
# dir.create(out, recursive=T)
# saveArchRProject(ArchRProj = proj, outputDirectory=out)

