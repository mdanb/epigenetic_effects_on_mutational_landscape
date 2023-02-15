library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
cores=args$cores
dataset = args$dataset
addArchRThreads(threads = cores)
addArchRGenome("hg19")

#h5disableFileLocking()
#h5enableFileLocking()
create_arrow_files <- function(fragment_paths) {
    sample_names = strsplit(fragment_paths, split="/")
    length = length(sample_names[[1]])
    sample_names = lapply(sample_names, "[", length)
    sample_names = sub(".fragments.*", "", sample_names)
    createArrowFiles(inputFiles = fragment_paths,
                     sampleNames = sample_names,
                     outputNames = sample_names,
	  	               minTSS = 4, 
                     minFrags = 1000,
                     addTileMat = T,
                     addGeneScoreMat = T,
                     force = F,
		                 cleanTmp = F,
		                 QCDir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/analysis/ArchR_proj")
				   
    #archp = ArchRProject(ArrowFiles = arrow_files, 
    #                     outputDirectory = output_dir,
    #                     copyArrows = F)
    #saveArchRProject(archp)
    #write.csv(archp@cellColMetadata, paste(output_dir, "metadata.csv", sep="/"))
 # }
}

root = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data"
if (dataset == "Bingren") {
  output_dir = paste(root, "arrow/Bingren", sep="/") 
  files_dir = paste(root, "bed_files/bingren_scATAC/migrated_to_hg19", sep="/")
} else if (dataset == "Shendure") {
  output_dir = paste(root, "arrow/Shendure", sep="/") 
  files_dir = paste(root, "bed_files/JShendure_scATAC", sep="/")
} else if (dataset == "Greenleaf_pbmc_bm") {
   output_dir = paste(root, "arrow/Greenleaf_pbmc_bm", sep="/")
   files_dir = paste(root, "bed_files/greenleaf_pbmc_bm_scATAC", 
                     sep="/")
} else if (dataset == "Greenleaf_brain") {
   output_dir = paste(root, "arrow/Greenleaf_brain", sep="/")
   files_dir = paste(root, "bed_files/greenleaf_brain_scATAC/migrated_to_hg19",
                     sep="/")
} else if (dataset == "Tsankov") {
   output_dir = paste(root, "arrow/Tsankov", sep="/")
   files_dir = paste(root, "bed_files/Tsankov_scATAC/migrated_to_hg19", sep="/")
} else if (dataset == "Yang_kidney") {
    output_dir = paste(root, "arrow/Yang_kidney", sep="/")
    files_dir = paste(root, "bed_files/yang_kidney_scATAC", sep="/")
}

setwd(root)
dir.create(root)
setwd(output_dir)

create_arrow_files(list.files(files_dir, full.names=T, pattern = "bgz$"))
