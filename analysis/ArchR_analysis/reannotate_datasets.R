library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--cores", type="integer"),
  make_option("--dataset", type="character"),
  make_option("--cluster", action="store_true"),
  make_option("--marker_genes", type="character", default=NULL)
)

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores
marker_genes = unlist(strsplit(marker_genes, split=","))

addArchRThreads(threads = cores)
addArchRGenome("hg19")

arrow_files_path = paste("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/arrow",
                          dataset, sep="/")

arrow_files = list.files(arrow_files_path, full.names=T, pattern="arrow")

dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj"
ArchR_proj <- loadArchRProject(dir)

cell_col_data = getCellColData(ArchR_proj)
dataset_subsetted_idx = cell_col_data[["dataset_per_cell"]] == dataset
ArchR_proj_subset = ArchR_proj[dataset_subsetted_idx]
cell_col_data_subset = getCellColData(ArchR_proj_subset)

ArchR_proj_subset <- addGeneScoreMatrix(ArchR_proj_subset)
ArchR_proj_subset <- addTileMatrix(ArchR_proj_subset)

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

ArchR_proj_subset <- addImputeWeights(ArchR_proj_subset, force=T)


#saveArchRProject(ArchRProj = ArchR_proj)

p <- plotEmbedding(
  ArchRProj = ArchR_proj, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95)
)

if (!("Clusters" %in% colnames(cell_col_data))) {
  ArchR_proj <- addClusters(
                                    input = ArchR_proj,
                                    reducedDims = "IterativeLSI",
                                    method = "Seurat",
                                    name = "Clusters",
                                    resolution = 0.3
                                  )
  saveArchRProject(ArchRProj = ArchR_proj)
}

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


pdf("../../figures/gene_marker_UMAPs_Tsankov_lung.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
dev.off()

p <- plotEmbedding(
  ArchRProj = ArchR_proj, 
  colorBy = "cellColData", 
  name = "Clusters", 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95)
)

plotPDF(p, name="clusters_UMAPs_Tsankov_lung.pdf", 
        ArchRProj = ArchR_proj, addDOC = FALSE)

