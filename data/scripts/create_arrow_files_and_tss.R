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

#h5disableFileLocking()
#h5enableFileLocking()
create_arrow_files_and_tss <- function(genome_version, fragment_paths, 
                                       output_dir) {
#if (!file.exists(paste(output_dir, "Save-ArchR-Project.rds", sep="/"))) {
    addArchRGenome(genome_version)
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

if (dataset == "bing_ren") {
  output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/bingren" 
  files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/bingren_scATAC/migrated_to_hg19"
  create_arrow_files_and_tss("hg19", list.files(output_dir, full.names=T,
                                                pattern = "bgz$"),
						output_dir)
						
} else if (dataset == "shendure") {
  output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/shendure"
  files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/JShendure_scATAC"
  create_arrow_files_and_tss("hg19", list.files(files_dir, full.names=T,
                                     pattern = "bgz$"), output_dir)
} else if (dataset == "greenleaf_pbmc_bm") {
   output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/greenleaf_pbmc_bm"
   files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_pbmc_bm_scATAC/migrated_to_hg19"
   create_arrow_files_and_tss("hg19", list.files(files_dir, full.names=T,
                               pattern = "bgz$"), output_dir)
} else if (dataset == "greenleaf_brain") {
   output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/greenleaf_brain"
   files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/greenleaf_brain_scATAC/migrated_to_hg19"
   create_arrow_files_and_tss("hg19", list.files(files_dir, full.names=T,
	                          pattern = "bgz$"), output_dir)

} else if (dataset == "tsankov") {
   output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/tsankov"
   files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/Tsankov_scATAC/migrated_to_hg19"
   create_arrow_files_and_tss("hg19", list.files(files_dir, full.names=T,
	                                   pattern = "bgz$"), output_dir)
} else if (dataset == "yang_kidney") {
    output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/arrow/yang_kidney"
    files_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/yang_kidney_scATAC"
    create_arrow_files_and_tss("hg19", list.files(files_dir, full.names=T,
	                          pattern = "bgz$"), output_dir)
}
