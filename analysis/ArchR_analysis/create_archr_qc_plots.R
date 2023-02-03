library(ArchR)

addArchRThreads(threads = 8)
addArchRGenome("hg19")

ArrowFiles = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow", 
                    			recursive=T,
                    			full.names=T,
                    			pattern=".arrow")
ArchR_proj <- ArchRProject(ArrowFiles = ArrowFiles,
			                     copyArrows=F,
    			      outputDirectory = "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/analysis/ArchR_analysis/ArchR_project_dir")

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
                            data=dataset_per_cell, name="dataset",
			                      cells=rownames(getCellColData(ArchRProj = ArchR_proj)))

ggplot(ArchR_proj$cellColData, aes(x = Sample, y = TSSEnrichment)) +
        geom_boxplot(aes(fill = dataset_per_cell), outlier.size=.2) +
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        xlab('Sample') + 
        ylab('TSSEnrichment')  

ggsave("../../figures/archr_tss_vs_sample.png")

ggplot(ArchR_proj$cellColData, aes(x = Sample, y = TSSEnrichment)) +
  geom_boxplot(aes(fill = dataset), outlier.size=.2) +
  theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Sample') + 
  ylab('Log(num fragments)')

ggsave("../../figures/archr_log_frags_vs_sample.png")
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