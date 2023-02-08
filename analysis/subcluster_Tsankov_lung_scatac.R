library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores

addArchRThreads(threads = cores)
addArchRGenome("hg19")

marker_genes = c('WT1','ITLN1','COL1A1','PECAM1','LYZ',
                 'CD3D','MSLN','KRT18','KRT5', 'VIM')

arrow_files = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/arrow/Tsankov", 
                         full.names=T, pattern="arrow")

dir = "archr_subcluster_Tsankov_lung"
dir.create(dir)
if (!file.exists("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/sub_cluster_Tsankov")) {
  sub_cluster_Tsankov <- ArchRProject(ArrowFiles = arrow_files, 
                                      outputDirectory = dir,
                                      copyArrows = FALSE)
  
  sub_cluster_Tsankov <- addGeneScoreMatrix(sub_cluster_Tsankov)
  sub_cluster_Tsankov <- addTileMatrix(sub_cluster_Tsankov)
  
  sub_cluster_Tsankov <- addIterativeLSI(
    ArchRProj = sub_cluster_Tsankov,
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
  
  sub_cluster_Tsankov <- addUMAP(ArchRProj = sub_cluster_Tsankov, 
                                 reducedDims = "IterativeLSI", 
                                 name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                                 metric = "cosine")
  saveArchRProject(ArchRProj = sub_cluster_Tsankov, 
                   outputDirectory = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/sub_cluster_Tsankov")
} else {
  sub_cluster_Tsankov = loadArchRProject("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/sub_cluster_Tsankov")
}

p <- plotEmbedding(
  ArchRProj = sub_cluster_Tsankov, 
  colorBy = "GeneScoreMatrix", 
  name = marker_genes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
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

pdf("../figures/gene_marker_UMAPs_Tsankov_lung.pdf")
cp = do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
dev.off()