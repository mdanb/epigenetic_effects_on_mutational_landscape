library(ArchR)
library(optparse)

option_list <- list( 
  make_option("--dataset", type="char")
)

args = parse_args(OptionParser(option_list=option_list))
dataset = args$dataset

addArchRThreads(threads = 8)

create_arrow_files_and_tss <- function(genome_version, fragment_paths, 
                                       output_dir) {
  if (!file.exists(paste(output_dir, "Save-ArchR-Project.rds", sep="/"))) {
    addArchRGenome(genome_version)
    sample_names = sub(".fragments.*", "", fragment_paths)
    arrow_files = createArrowFiles(inputFiles = fragment_paths,
                                   sampleNames = sample_names,
                                   minTSS = 4, 
                                   minFrags = 1000, 
                                   addTileMat = F,
                                   addGeneScoreMat = F,
                                   force = F)
    archp = ArchRProject(ArrowFiles = arrow_files, 
                         outputDirectory = output_dir,
                         copyArrows = FALSE)
    saveArchRProject(archp)
    write.csv(archp@cellColMetadata, paste(output_dir, "metadata.csv", sep="/"))
  }
}

if (dataset == "bing_ren") {
  output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files/JShendure_scATAC"
  create_arrow_files_and_tss("hg38", list.files(output_dir), output_dir)
} else if (dataset == "shendure") {
  output_dir = "/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/bed_files"
  create_arrow_files_and_tss("hg19", list.files(output_dir), output_dir)
}
