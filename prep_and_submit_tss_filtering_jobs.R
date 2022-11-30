library(optparse)
library(tidyverse)

dir.create("tsse_filtered_run_commands/for_ML_input")
dir.create("tsse_filtered_run_commands/for_ML_input/log_dir")

option_list <- list( 
  make_option("--bing_ren", action="store_true", default = FALSE),
  make_option("--shendure", action="store_true", default = FALSE),
  make_option("--tsankov", action="store_true", default = FALSE),
  make_option("--fragment_count_range", type="character")
)

args = parse_args(OptionParser(option_list=option_list))

# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--fragment_count_range=100000,150000,200000",
#                      "--bing_ren",
#                      "--shendure",
#                      "--tsankov"))


bing_ren = args$bing_ren
shendure = args$shendure
tsankov = args$tsankov
fragment_count_range = args$fragment_count_range

create_and_submit_job_scipts <- function(bing_ren, shendure, tsankov,
                                         fragment_count_range) {
  metadatas = list()
  raw_data_dirs = c()
  
  if (bing_ren) {
    metadata_bingren = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", 
                                  header=TRUE)
    colnames(metadata_bingren)[grep("cell.type", colnames(metadata_bingren))] = "cell_type"
    metadata_bingren = metadata_bingren %>% filter(grepl("SM",tissue))
    metadata_bingren["files_pattern"] = str_extract(metadata_bingren[["tissue"]], 
                                                    ".*_SM")
    metadatas = c(metadatas, list(metadata_bingren))
    
    bingren_raw_dir = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                            "raw_dir", "bed_files", sep="/")
    raw_data_dirs = append(raw_data_dirs, bingren_raw_dir)
  }
  
  if (shendure) {
    metadata_shendure = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                                   sep="\t",
                                   header=TRUE)
    colnames(metadata_shendure)[grep("tissue", colnames(metadata_shendure))] = "files_pattern"
    metadatas = c(metadatas, list(metadata_shendure))
    
    shendure_raw_dir = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                             "raw_dir", "bed_files", "JShendure_scATAC", sep="/")
    raw_data_dirs = append(raw_data_dirs, shendure_raw_dir)
  }
  
  if (tsankov) {
    metadata_tsankov_proximal = 
      read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/tsankov_lung_proximal_barcode_annotation.csv")
    metadata_tsankov_distal = 
      read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv")
    metadata_tsankov_distal["files_pattern"] = "RPL"
    metadata_tsankov_proximal["files_pattern"] = "IC"
    colnames(metadata_tsankov_proximal)[grep("celltypes", 
                                             colnames(metadata_tsankov_proximal))] = "cell_type"
    colnames(metadata_tsankov_distal)[grep("celltypes", 
                                           colnames(metadata_tsankov_distal))] = "cell_type"
    metadatas = c(metadatas, list(metadata_tsankov_proximal))
    metadatas = c(metadatas, list(metadata_tsankov_distal))
    
    tsankov_raw_dir = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                            "raw_dir", "bed_files", "Tsankov_scATAC", sep="/")
    raw_data_dirs = append(raw_data_dirs, tsankov_raw_dir)
  }
  
  idx = 1
  for (metadata in metadatas) {
    metadata = as.data.frame(metadata)
    raw_data_dir = raw_data_dirs[idx]
    mem = 32
    if (grepl("BingRen", raw_data_dir)) {
      dataset = "bing_ren"
    } else if ("Shendure" %in% raw_data_dir) {
      dataset = "shendure"
      mem = 64
    } else {
      dataset = "tsankov"
    }
    for (pattern in unique(metadata[["files_pattern"]])) {
      print(paste("creating and submitting job for", pattern, sep=" "))
      files = setdiff(list.files(raw_data_dir, pattern=pattern), 
                      list.dirs(raw_data_dir, 
                                recursive = FALSE, 
                                full.names = FALSE))
      num_files = length(files)
      cell_types = metadata %>% 
                   filter(files_pattern == pattern) %>%
                   select(cell_type)
                      
      cell_types = unique(cell_types)
      job_script_filename = paste(paste(pattern, dataset, sep="_"), "sh", 
                                  sep=".")
      job_script_path = paste("tsse_filtered_run_commands/for_ML_input", 
                              job_script_filename, sep="/")
      file_conn = file(paste(job_script_path, sep="/"))
      mem_req = paste0("#$ -l h_vmem=", mem * num_files, "G")
      
      pattern_dataset = paste(pattern, dataset, sep="_")
      job_name = paste("#$ -N", pattern_dataset)
      output_log_filename = paste("#$ -o /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin",
                                  "tsse_filtered_run_commands", "for_ML_input", 
                                  "log_dir", 
                                  paste(pattern_dataset, "out", sep="."), 
                                  sep = "/")
      error_log_filename = paste("#$ -e /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin", 
                                 "tsse_filtered_run_commands", "for_ML_input", 
                                 "log_dir",
                                 paste(pattern_dataset, "err", sep="."), 
                                 sep = "/")
      
      cell_types = paste(as.list(cell_types)$cell_type, collapse=";")
      
      command_args = paste0("--cell_types='", cell_types, "' ",
                            "--dataset='", dataset, "' ",
                            "--top_tsse_fragment_count_range='", fragment_count_range, "' ",
                            "--files_pattern='", pattern, "' ",
                            "--cores=", 1)
      
      rscript_command = paste("Rscript /ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/create_tsse_filtered_count_overlaps.R",
                              command_args)
      
      writeLines(c("#!/bin/bash", 
                   "#$ -cwd", 
                   "#$ -l h_rt=24:00:00", 
                   "#$ -l os=RedHat7", 
                   "#$ -pe smp 1", 
                   mem_req, job_name, 
                   output_log_filename, 
                   error_log_filename, 
                   "#$ -binding linear:1",
                   "#$ -M bgiotti@broadinstitute.org",
                   "#$ -m bea",
                   "",
                   "source /broad/software/scripts/useuse",
                   "use .anaconda3-5.3.1",
                   "source activate /home/unix/bgiotti/conda/coo",
                   "",
                   rscript_command), 
                file_conn)
      system(paste("qsub", job_script_path))
    }
    idx = idx + 1
  }
}

create_and_submit_job_scipts(bing_ren, shendure, tsankov, fragment_count_range)
