library(ArchR)
library(optparse)

option_list <- list(
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores

addArchRThreads(threads = cores)
addArchRGenome("hg19")

ArrowFiles = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj/ArrowFiles", 
                        recursive=T,
                        full.names=T,
                        pattern=".arrow")
ArchR_proj <- ArchRProject(ArrowFiles = ArrowFiles,
                           copyArrows=F,
                           outputDirectory = 
                          "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj")

dataset = unlist(lapply(strsplit(ArrowFiles, split = "/"), "[", 8))
sample_col_data = getSampleColData(ArchRProj =  ArchR_proj)
ArchR_proj = addSampleColData(ArchRProj = ArchR_proj, 
                              data=dataset, name="dataset", 
                              samples=rownames(sample_col_data))

samples = rownames(sample_col_data)
cell_col_data = getCellColData(ArchRProj = ArchR_proj)
samples_per_cell = unlist(lapply(strsplit(rownames(cell_col_data), 
                                          split = "#"), "[", 1))
idx = match(samples_per_cell, samples)
dataset_per_cell = dataset[idx]

ArchR_proj = addCellColData(ArchRProj = ArchR_proj,
                            data=dataset_per_cell, name="dataset_per_cell",
                            cells=rownames(cell_col_data))

# ArchR_proj = addCellColData(ArchRProj = ArchR_proj,
#                             data=metadata_per_cell, name="metadata_per_cell",
#                             cells=rownames(cell_col_data))
saveArchRProject(ArchRProj = ArchR_proj)
