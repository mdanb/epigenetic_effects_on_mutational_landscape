library(ArchR)
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(dplyr)

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character", default="all"),
  make_option("--metadata_for_celltype_fn", type="character"),
  make_option("--sep_for_metadata", type="character", default=","),
  make_option("--cell_type_col_in_metadata", type="character"),
  # make_option("--cell_name_col_in_metadata", type="character"),
  make_option("--cluster", action="store_true", default=F),
  make_option("--cluster_res", type="double", default=1.2),
  make_option("--plot_cell_types", action="store_true", default=F),
  make_option("--tissue", type="character", default="all"),
  make_option("--nfrags_filter", type="integer", default=1),
  make_option("--tss_filter", type="integer", default=0),
  make_option("--tss_percentile", type="double", default=NULL),
  make_option("--nfrags_percentile", type="double", default=NULL),
  make_option("--filter_per_cell_type", action="store_true", default=FALSE),
  make_option("--cell_types", type="character", default="all"),
  make_option("--marker_genes", type="character", default=NULL),
  make_option("--min_cells_per_cell_type", type="integer", default=0)
)

# args = parse_args(OptionParser(option_list=option_list), args=
#                     c("--cores=8",
#                       "--dataset=Tsankov",
#                       "--metadata_for_celltype_fn=combined_distal_proximal.csv",
#                       "--sep_for_metadata=,",
#                       "--filter_per_cell_type",
#                       "--cell_type_col_in_metadata=celltypes",
#                       "--cluster",
#                       "--cluster_res=1",
#                       "--tissue=RPL",
#                       "--nfrags_filter=1",
#                       "--tss_filter=0",
#                       "--cell_types=all",
#                       "--min_cells_per_cell_type=1")
# )

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

print("Collecting cmd line args")
# args = parse_args(OptionParser(option_list=option_list))
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
  print("Done loading existing project")
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

if (cluster) {
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

  fn = paste("clusters_UMAPs", setting, sep="_")
  fn = paste0(fn, ".pdf")
  print(paste("saving", fn))
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}

if (plot_cell_types) {
  p <- plotEmbedding(
        ArchRProj = proj, 
        colorBy = "cellColData", 
        name = "cell_type", 
        embedding = "UMAP",
        quantCut = c(0.01, 0.95))

  fn = paste("cell_type_UMAP", setting, sep="_")
  fn = paste0(fn, ".pdf")
  plotPDF(p, name=fn, ArchRProj = proj, addDOC = FALSE)
}

if (!(is.null(marker_genes))) {
  proj <- addImputeWeights(proj)
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
    baseSize=20
  )
  if (length(marker_genes) == 1) {
    p = list(p)
  }
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 3) +
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
  pdf(fp, width = 50, height = 50)
  do.call(cowplot::plot_grid, c(list(ncol = 5), p2))
  dev.off()
}
########## Post visual inspection ##########
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
}

saveArchRProject(ArchRProj = proj)



# out = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis", 
#             "ArchR_per_dataset", dataset, tissue, sep="/")
# dir.create(out, recursive=T)
# saveArchRProject(ArchRProj = proj, outputDirectory=out)

