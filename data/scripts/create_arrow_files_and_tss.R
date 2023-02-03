library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--dataset", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--dataset=bing_ren"))
dataset = args$dataset
# dataset="bing_ren"
#Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE)
#Sys.setenv(RHDF5_USE_FILE_LOCKING=FALSE)
addArchRThreads(threads = 8)
addArchRGenome("hg19")

#h5disableFileLocking()
#h5enableFileLocking()
create_arrow_files <- function(fragment_paths, 
                               output_dir) {
#if (!file.exists(paste(output_dir, "Save-ArchR-Project.rds", sep="/"))) {
    sample_names = strsplit(fragment_paths, split="/")
    length = length(sample_names[[1]])
    sample_names = lapply(sample_names, "[", length)
    sample_names = sub(".fragments.*", "", sample_names)
    createArrowFiles(inputFiles = fragment_paths,
                     sampleNames = sample_names,
                     outputNames = paste(output_dir, sample_names, sep="/"),
	  	               minTSS = 4, 
                     minFrags = 1000,
                     addTileMat = F,
                     addGeneScoreMat = F,
                     force = F,
		     QCDir = "../../analysis/ArchR_analysis/QualityControl")
				   
    #archp = ArchRProject(ArrowFiles = arrow_files, 
    #                     outputDirectory = output_dir,
    #                     copyArrows = F)
    #saveArchRProject(archp)
    #write.csv(archp@cellColMetadata, paste(output_dir, "metadata.csv", sep="/"))
 # }
}

root = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data"
if (dataset == "Bingren") {
  output_dir = paste(root, "arrow/bingren", sep="/") 
  files_dir = paste(root, "bed_files/bingren_scATAC/migrated_to_hg19", sep="/")
} else if (dataset == "Shendure") {
  output_dir = paste(root, "arrow/shendure", sep="/") 
  files_dir = paste(root, "bed_files/JShendure_scATAC", sep="/")
} else if (dataset == "Greenleaf_pbmc_bm") {
   output_dir = paste(root, "arrow/greenleaf_pbmc_bm", sep="/")
   files_dir = paste(root, "bed_files/greenleaf_pbmc_bm_scATAC/migrated_to_hg19", 
                     sep="/")
} else if (dataset == "Greenleaf_brain") {
   output_dir = paste(root, "arrow/greenleaf_brain", sep="/")
   files_dir = paste(root, "bed_files/greenleaf_brain_scATAC/migrated_to_hg19",
                     sep="/")
} else if (dataset == "Tsankov") {
   output_dir = paste(root, "arrow/tsankov", sep="/")
   files_dir = paste(root, "bed_files/Tsankov_scATAC/migrated_to_hg19", sep="/")
} else if (dataset == "Yang_kidney") {
    output_dir = paste(root, "arrow/yang_kidney", sep="/")
    files_dir = paste(root, "bed_files/yang_kidney_scATAC", sep="/")
}

create_arrow_files(list.files(files_dir, full.names=T, pattern = "bgz$"),  
                   output_dir)