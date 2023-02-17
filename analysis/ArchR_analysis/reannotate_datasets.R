library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character"),
  make_option("--cluster", action="store_true"),
  make_option("--cluster_res", type="integer"),
  make_option("--marker_genes", type="character", default=NULL)
)

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores
dataset = args$dataset
cluster = args$cluster
cluster_res = args$cluster_res
marker_genes = args$marker_genes

addArchRThreads(threads = cores)
addArchRGenome("hg19")

# arrow_files_path = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/arrow",
#                           dataset, sep="/")

# arrow_files = list.files(arrow_files_path, full.names=T, pattern="arrow")
# markers = c('WT1','ITLN1','COL1A1','PECAM1','LYZ','CD3D','MSLN','KRT18','KRT5','VIM')
dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj"
ArchR_proj <- loadArchRProject(dir)

cell_col_data = getCellColData(ArchR_proj)
dataset_subsetted_idx = cell_col_data[["dataset_per_cell"]] == dataset
ArchR_proj_subset = ArchR_proj[dataset_subsetted_idx]
cell_col_data_subset = getCellColData(ArchR_proj_subset)

# tryCatch({
# getMatrixFromProject(ArchR_proj_subset, useMatrix = "GeneScoreMatrix")
# }, error = function(err) {
ArchR_proj_subset <- addGeneScoreMatrix(ArchR_proj_subset, force=T)
# })

# tryCatch({
  # getMatrixFromProject(ArchR_proj_subset, useMatrix = "TileMatrix")
# }, error = function(err) {
ArchR_proj_subset <- addTileMatrix(ArchR_proj_subset, force=T)
# })


ArchR_proj_subset <- addIterativeLSI(
  ArchRProj = ArchR_proj_subset,
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
  force=T
)

ArchR_proj_subset <- addUMAP(ArchRProj = ArchR_proj_subset, 
                             reducedDims = "IterativeLSI", 
                             name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                             metric = "cosine",
                             force=T)

ArchR_proj_subset <- addImputeWeights(ArchR_proj_subset)


#saveArchRProject(ArchRProj = ArchR_proj)
if (!(is.null(marker_genes))) {
  marker_genes = unlist(strsplit(marker_genes, split=","))
  p <- plotEmbedding(
    ArchRProj = ArchR_proj_subset, 
    colorBy = "GeneScoreMatrix", 
    name = marker_genes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95)
  )
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  path = "../../figures"
  fn = paste0("gene_marker_UMAPs_", dataset, ".pdf")
  fp = paste(path, fn, sep="/")
  pdf(fp)
  do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
  dev.off()
}

if (cluster) {
# if (!("Clusters" %in% colnames(cell_col_data))) {
  ArchR_proj_subset <- addClusters(
                                    input = ArchR_proj_subset,
                                    reducedDims = "IterativeLSI",
                                    method = "Seurat",
                                    name = "Clusters",
                                    resolution = cluster_res,
                                    force=T)
  p <- plotEmbedding(
    ArchRProj = ArchR_proj_subset, 
    colorBy = "cellColData", 
    name = "Clusters", 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95)
  )
  
  fn = paste0("clusters_UMAPs_", dataset, ".pdf")
  plotPDF(p, name=fn, ArchRProj = ArchR_proj_subset, addDOC = FALSE)
}

out = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis", 
            "ArchR_per_dataset", dataset, sep="/")
dir.create(out, recursive=T)
saveArchRProject(ArchRProj = ArchR_proj_subset, outputDirectory=out)

