library(ArchR)

addArchRThreads(threads = 8)
addArchRGenome("hg19")

ArrowFiles = list.files(
ArchR_proj <- ArchRProject(ArrowFiles = ArrowFiles, 
    			   outputDirectory = "ArchR_analysis",
     			   copyArrows = TRUE)


