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
		                 cleanTmp = T,
		               QCDir = "../../../analysis/ArchR_analysis/ArchR_projects/QC")
				   
    #archp = ArchRProject(ArrowFiles = arrow_files, 
    #                     outputDirectory = output_dir,
    #                     copyArrows = F)
    #saveArchRProject(archp)
    #write.csv(archp@cellColMetadata, paste(output_dir, "metadata.csv", sep="/"))
 # }
}

root = ".."
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
} else if (dataset == "Greenleaf_colon") {
  output_dir = paste(root, "arrow/Greenleaf_colon", sep="/")
  files_dir = paste(root, "bed_files/greenleaf_colon_scATAC/migrated_to_hg19", 
                    sep="/")
} else if (dataset == "Rawlins_fetal_lung") {
  output_dir = paste(root, "arrow/Rawlins_fetal_lung", sep="/")
  files_dir = paste(root, 
                    "bed_files/rawlins_fetal_lung_scATAC/migrated_to_hg19", 
                    sep="/")
}

dir.create(output_dir)
setwd(output_dir)
files_dir = paste(files_dir, sep="/")
files = list.files(files_dir, full.names=T, pattern = "bgz$")
if (dataset == "Yang_kidney") {
  files = list.files(files_dir, full.names=T, pattern = "gz$")
}

create_arrow_files(files)
