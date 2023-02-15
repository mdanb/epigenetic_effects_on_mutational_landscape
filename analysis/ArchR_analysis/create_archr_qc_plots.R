library(ArchR)
library(optparse)
library(tibble)

option_list <- list(
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
cores = args$cores

addArchRThreads(threads = cores)
addArchRGenome("hg19")

ArrowFiles = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/arrow", 
                    			recursive=T,
                    			full.names=T,
                    			pattern=".arrow")
ArchR_proj <- ArchRProject(ArrowFiles = ArrowFiles,
			                     copyArrows=F,
    			      outputDirectory = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj")

dataset = unlist(lapply(strsplit(ArrowFiles, split = "/"), "[", 8))
ArchR_proj = addSampleColData(ArchRProj = ArchR_proj, 
			                        data=dataset, name="dataset", 
		    	                    samples=rownames(getSampleColData(ArchRProj = 
		    	                                                        ArchR_proj)))

sample_col_data = getSampleColData(ArchRProj = ArchR_proj)
#bingren = sample_col_data[sample_col_data$dataset == "bingren", ]
#samples = rownames(bingren)
samples = rownames(sample_col_data)
samples_per_cell = unlist(lapply(strsplit(rownames(getCellColData(ArchRProj = ArchR_proj)), 
			  split = "#"), "[", 1))
idx = match(samples_per_cell, samples)
dataset_per_cell = dataset[idx]
ArchR_proj = addCellColData(ArchRProj = ArchR_proj,
                            data=dataset_per_cell, name="dataset_per_cell",
			                      cells=rownames(getCellColData(ArchRProj = ArchR_proj)))

ggplot(as_tibble(getCellColData(ArchR_proj)), 
        aes(x = Sample, y = TSSEnrichment)) +
        geom_boxplot(aes(fill = dataset_per_cell), outlier.size=.2) +
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab('Sample') + 
        ylab('TSSEnrichment') +
        labs(fill="dataset")

ggsave("../../figures/archr_tss_vs_sample.png", width=20)

ggplot(as_tibble(getCellColData(ArchR_proj)), 
       aes(x = Sample, y = log10(nFrags))) +
       geom_boxplot(aes(fill = dataset_per_cell), outlier.size=.2) +
       theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       xlab('Sample') + 
       ylab('Log(num fragments)') +
       labs(fill="dataset")

ggsave("../../figures/archr_log_frags_vs_sample.png", width=20)
saveArchRProject(ArchRProj = ArchR_proj, load = FALSE)
#bingren_cells = ArchR_proj

# tss_p <-plotGroups(
#               ArchRProj = ArchR_proj, 
#               groupBy = "Sample", 
#               colorBy = cellColData,
#               name = "TSSEnrichment",
#               plotAs = "violin",
#               alpha = 0.4,
#               addBoxPlot = TRUE,
#               ratioYX=1/3)

# num_fragments_p <- plotGroups(
#                           ArchRProj = ArchR_proj, 
#                           groupBy = "Sample", 
#                           colorBy = cellColData,
#                           name = "num",
#                           plotAs = "violin",
#                           alpha = 0.4,
#                           addBoxPlot = TRUE,
#                           ratioYX=1/3)

# plotPDF(tss_p, name = "QC-Sample-Statistics.pdf", ArchRProj = ArchR_proj, 
#         addDOC = FALSE, width = 20, height = 10)
# plotPDF(num_fragments_p, name = "QC-Sample-Statistics.pdf", 
#         ArchRProj = ArchR_proj, 
#         addDOC = FALSE, width = 20, height = 10)
