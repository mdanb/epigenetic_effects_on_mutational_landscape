library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(parallel)
library(optparse)
source("create_count_overlaps_utils.R")

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--cell_number_filter", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
cell_number_filter = args$cell_number_filter
dataset = args$dataset
# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 1) {
#   stop("One argument must be supplied", call.=FALSE)
# } else {
#   CELL_NUMBER_FILTER = as.integer(args[1])
# } 

compute_count_overlaps <- function(sample, interval_ranges) {
  grl_in = sample %>%                                             
    group_split(cell_type) %>%                             
    lapply(makeGRangesFromDataFrame, keep.extra.columns=TRUE)
  grl = GRangesList(grl_in)
  count_overlaps = lapply(grl, function(x) countOverlaps (interval_ranges, x))
  cell_types = unlist(lapply(grl, function(x) unique(x$cell_type)))
  names(count_overlaps) = cell_types
  return(count_overlaps)
}

create_count_overlaps_file <- function(file, cell_number_filter, metadata,
                                       interval_ranges, chain) {
  filename = paste("count_filter", cell_number_filter,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  if (!file.exists(filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
    sample_name = get_sample_name(file)
    filtered_metadata = filter_metadata_by_sample_name(sample_name, metadata)
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(filtered_metadata,
                                                                  "cellID",
                                                                  "\\+")
    sample <- filter_sample_to_contain_only_cells_in_metadata(sample,
                                                    sample_barcodes_in_metadata)
    sample <- get_sample_cell_types(sample, sample_barcodes_in_metadata, 
                                    filtered_metadata)
    
    # if (!is.na(cell_types_for_tsse_subsampling)) {
    #   search_name = paste(str_to_title(sample_name), sample$cell_type)
    #   sample = sample[search_name %in% cell_types_for_tsse_subsampling, ]
    #   full_cell_type_name_in_metadata = paste(str_to_title(sample_name), 
    #                                           filtered_metadata$cell.type)
    #   filtered_metadata = filtered_metadata[full_cell_type_name_in_metadata %in%
    #                                           cell_types_for_tsse_subsampling, ]
    # }
    
    sample = as_tibble(migrate_bed_file_to_hg37(sample, chain))
    sample_cell_counts_filename = paste0("cell_counts_", file_path_sans_ext(
                                          file, TRUE), ".rds")
    sample_cell_counts_file_path = paste("processed_data/cell_counts_per_sample", 
                                         sample_cell_counts_filename, sep="/")
    
    # if (!is.na(cell_types_for_tsse_subsampling)) {
    #   add_cumulative_num_fragments_per_cell_to_metadata
    #   sample <- filter_sample_by_tss(sample, subsampling_frag_count, 
    #                                  filtered_metadata)
    #   counts_per_cell_type <- readRDS(sample_cell_counts_file_path)
    # }
    # else {
    #   counts_per_cell_type <- get_num_cells_per_sample(sample)
    #   saveRDS(counts_per_cell_type, sample_cell_counts_file_path)
    # }
    sample <- filter_sample_by_cell_number(sample,
                                           counts_per_cell_type, 
                                           cell_number_filter)
    
    count_overlaps <- compute_count_overlaps(sample, interval_ranges)
    saveRDS(count_overlaps, filepath)
  }
  
  
  
  # if (!is.na(top_tsse_fragment_count_subsampling_range)) {
  #   for (fragment_count in top_tsse_fragment_count_subsampling_range) {
  #     filename = paste("tsse_filter", fragment_count, filename, sep="_")
  #     filepath = paste("processed_data/count_overlap_data/tsse_filtered", 
  #                      filename, sep="/")
  #     create_count_overlaps_file_helper(file, filepath, chain, metadata, 
  #                                       cell_number_filter, interval_ranges,
  #                                       cell_types_for_tsse_subsampling,
  #                                       fragment_count)
  #   }
  # }
  # else {
  #   filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  #   create_count_overlaps_file_helper(file, filepath, chain, metadata, 
  #                                     cell_number_filter, interval_ranges)
  # }
  
  # # subdivided_filename = paste("count_filter", CELL_NUMBER_FILTER,  
  # #                             "subdivided_count_overlaps", 
  # #                             paste(file_path_sans_ext(file, TRUE),
  # #                                   "rds", sep="."), sep="_")
  # # subdivided_filepath = paste("processed_data/count_overlap_data", 
  # #                             subdivided_filename, 
  # #                             sep="/")
  # # if (!file.exists(filepath) || !file.exists(subdivided_filepath)) {
  # if (!file.exists(filepath)) {
  #   print(paste("Processing", file, sep= " "))
  #   sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
  #   sample = migrate_bed_file_to_hg37(sample, chain)
  #   sample_name = get_sample_name(file)
  #   sample <- get_sample_cell_types(sample, sample_name, metadata)
  #   counts_per_cell_type <- get_and_save_num_cells_per_sample(sample, file)
  #   sample <- filter_sample_by_cell_number(sample,
  #                                          counts_per_cell_type, 
  #                                          cell_number_filter)
  #   # if (!file.exists(filename)) {
  #   count_overlaps <- compute_count_overlaps(sample, interval_ranges)
  #   saveRDS(count_overlaps, filepath)
  #   # }
  #   # if (!file.exists(subdivided_filename)) {
  #   #   subdivided_interval_ranges = subdivideGRanges(interval.ranges, 
  #   #                                                 subsize = 200000)
  #   #     # names(subdivided_interval_ranges) = rep(colnames(count_overlaps), each=5)
  #   #   
  #   #   count_overlaps_subdivided <- compute_count_overlaps(sample, 
  #   #                                                  subdivided_interval_ranges)
  #   #   # colnames(count_overlaps) = rep(colnames(count_overlaps), each=5)
  #   #   
  #   #   saveRDS(count_overlaps_subdivided, subdivided_filepath)
  #   # }
}

# create_count_overlaps_file_helper <- function(file, filepath, chain, metadata,
#                                               cell_number_filter, 
#                                               interval_ranges) {
#   if (!file.exists(filepath)) {
#     print(paste("Processing", file, sep= " "))
#     sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
#     sample_name = get_sample_name(file)
#     
#     filtered_metadata = metadata[metadata$sample %like% sample_name, ]
#     sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(filtered_metadata,
#                                                                   "cellID",
#                                                                   "\\+")
#     sample <- filter_sample_to_contain_only_cells_in_metadata(sample,
#                                                               sample_barcodes_in_metadata)
#     sample <- get_sample_cell_types(sample, filtered_metadata)
#     
#     if (!is.na(cell_types_for_tsse_subsampling)) {
#       search_name = paste(str_to_title(sample_name), sample$cell_type)
#       sample = sample[search_name %in% cell_types_for_tsse_subsampling, ]
#       full_cell_type_name_in_metadata = paste(str_to_title(sample_name), 
#                                               filtered_metadata$cell.type)
#       filtered_metadata = filtered_metadata[full_cell_type_name_in_metadata %in%
#                                             cell_types_for_tsse_subsampling, ]
#     }
#     
#     sample = as_tibble(migrate_bed_file_to_hg37(sample, chain))
#     sample_cell_counts_filename = paste0("cell_counts_", file_path_sans_ext(
#                                          file, TRUE), 
#                                         ".rds")
#     sample_cell_counts_file_path = paste("processed_data/cell_counts_per_sample", 
#                                          sample_cell_counts_filename, sep="/")
#     if (!is.na(cell_types_for_tsse_subsampling)) {
#       add_cumulative_num_fragments_per_cell_to_metadata
#       sample <- filter_sample_by_tss(sample, subsampling_frag_count, 
#                                      filtered_metadata)
#       counts_per_cell_type <- readRDS(sample_cell_counts_file_path)
#     }
#     else {
#       counts_per_cell_type <- get_num_cells_per_sample(sample)
#       saveRDS(counts_per_cell_type, sample_cell_counts_file_path)
#     }
#     sample <- filter_sample_by_cell_number(sample,
#                                            counts_per_cell_type, 
#                                            cell_number_filter)
#     
#     count_overlaps <- compute_count_overlaps(sample, interval_ranges)
#     saveRDS(count_overlaps, filepath)
#   }
# }

create_count_overlaps_file_shendure <- function(file, cell_number_filter, 
                                                metadata, interval_ranges) {
  filename = paste("Shendure_count_filter", cell_number_filter,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
                                           "rds", sep="."), sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  if (!file.exists(filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste("raw_dir", "bed_files", "JShendure_scATAC", file, 
                          sep="/"), format="bed")
    sample_name = get_sample_name_shendure(file)
    filtered_metadata = filter_metadata_by_sample_name(sample_name, metadata)
    sample_barcodes_in_metadata = filtered_metadata[["cell"]]
    sample <- filter_sample_to_contain_only_cells_in_metadata(sample,
                                                              sample_barcodes_in_metadata)
    sample <- get_sample_cell_types_shendure(sample, filtered_metadata)
    counts_per_cell_type <- get_and_save_num_cells_per_sample(sample, file)
    sample <- filter_sample_by_cell_number(sample,
                                           counts_per_cell_type, 
                                           cell_number_filter)
    if (!file.exists(filename)) {
      count_overlaps <- compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
  }
}

create_count_overlaps_file_tsankov <- function(file, cell_number_filter, 
                                               metadata, interval_ranges, 
                                               chain) {
  filename = paste(file_path_sans_ext(file, TRUE), "rds", sep=".")
  filename = unlist(strsplit(filename, split = "_"))
  filename = paste(filename[2:length(filename)], collapse="_")
  filename = paste("Tsankov_count_filter", cell_number_filter,  
                   "count_overlaps", filename, sep="_")
  filepath = paste("processed_data/count_overlap_data", filename, sep="/")
  if (!file.exists(filepath)) {
    print(paste("Processing", file, sep= " "))
    sample = import(paste0(paste("raw_dir", "bed_files", "Tsankov_scATAC", 
                                 substr(file, 1, nchar(file) - 3), sep="/"), 
                           "tsv"),
                    format="bed")
    sample = migrate_bed_file_to_hg37(sample, chain)
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(metadata, 
                                                                  "X",
                                                                  "#")
    sample <- filter_sample_to_contain_only_cells_in_metadata(sample,
                                                              sample_barcodes_in_metadata)
    sample <- get_sample_cell_types_tsankov(sample, metadata)
    counts_per_cell_type <- get_and_save_num_cells_per_sample(sample, file)
    sample <- filter_sample_by_cell_number(sample,
                                           counts_per_cell_type, 
                                           cell_number_filter)
    if (!file.exists(filename)) {
      count_overlaps = compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
  }
}

# create_count_overlaps_file_per_tsse_count <- function(count, 
#                                           metadata_with_fragment_counts_per_cell_type,
#                                           migrated_fragments,
#                                           sample_names, interval_ranges,
#                                           cell_types_to_consider, 
#                                           cell_number_filter) {
#   filename = paste("count_filter", cell_number_filter,  
#                    "count_overlaps", "frag_count_filter", count, sep="_")
#   filename = paste(filename, "rds", sep=".")
#   tissue_name = get_tissue_name(files[1])
#   filepath = paste("processed_data/count_overlap_data/tsse_filtered", 
#                    tissue_name, sep="/")   
#   dir.create(filepath)
#   filepath = paste(filepath, filename, sep="/")
#   if (!file.exists(filepath)) {
#     debug_message = paste0("Creating count overlaps for num fragments = ", 
#                            count)
#     print(debug_message)
#     count_filtered_metadata_with_fragment_counts_per_cell_type = 
#       mclapply(metadata_with_fragment_counts_per_cell_type, 
#                filter_metadata_by_fragment_count, count, 
#                mc.cores=8)
#     # top_barcodes = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type,
#     #                       function(x) x[["cell_barcode"]])
#     fragments_from_top_cells = mclapply(seq_along(count_filtered_metadata_with_fragment_counts_per_cell_type),
#                                         get_fragments_from_top_cells,
#                                         count_filtered_metadata_with_fragment_counts_per_cell_type,
#                                         migrated_fragments,
#                                         sample_names, mc.cores=8)
#     fragments_from_top_cells = mclapply(fragments_from_top_cells, sample, count,
#                                         mc.cores=8)
#     count_overlaps = mclapply(fragments_from_top_cells, 
#                               function(x) countOverlaps(interval_ranges, x),
#                               mc.cores=8) 
#     count_overlaps = as_tibble(count_overlaps, .name_repair="universal")
#     colnames(count_overlaps) = cell_types_to_consider
#     saveRDS(count_overlaps, filepath)
#   }
# }

filter_sample_by_cell_number <- function(sample, 
                                         counts_per_cell_type, 
                                         cell_number_filter) {
  cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= 
                                          cell_number_filter, ]$cell_type
  sample = sample %>%                                             
           filter(cell_type %in% cell_type_keep)
  return(sample)
}

filter_sample_to_contain_only_cells_in_metadata <- function(sample,
                                                            sample_barcodes_in_metadata) {
  sample = sample[sample$name %in% sample_barcodes_in_metadata]
  return(sample)
}

get_num_cells_per_sample <- function(sample) {
  counts_per_cell_type = sample %>%                               
                         group_by(cell_type) %>%                 
                         summarize(n_cells = n_distinct(name))
  return(counts_per_cell_type)
}

get_sample_cell_types <- function(fragments, sample_barcodes_in_metadata,
                                  filtered_metadata) {
  sample_idx_in_metadata = match(fragments$name, 
                                 sample_barcodes_in_metadata)
  fragments$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                           "cell.type"]
  return(fragments)
}

get_sample_cell_types_shendure <- function(sample, sample_name, 
                                           filtered_metadata) {
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                       "cell_type"]
  return(as.tibble(sample))
}

get_sample_cell_types_tsankov <- function(sample, metadata) {
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = metadata[unlist(sample_idx_in_metadata),
                              "celltypes"]
  return(as.tibble(sample))
}

# metadata_tsankov_proximal = read.table("raw_dir/metadata/tsankov_lung_proximal_barcode_annotation.csv",
#                                        sep=",",
#                                        header=TRUE)
# metadata_tsankov_distal = read.table("raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv",
#                                        sep=",",
#                                        header=TRUE)
# colnames(metadata_tsankov_proximal)[2] <- "celltypes"
# colnames(metadata_tsankov_distal)[2] <- "celltypes"

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

dir.create("processed_data/count_overlap_data", recursive=TRUE)                                       
dir.create("processed_data/cell_counts_per_sample")                                       

if (dataset == "bing_ren") {
  metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t",
                        header=T)
  files = setdiff(list.files("raw_dir/bed_files/"), 
                  list.dirs("raw_dir/bed_files", recursive = FALSE, 
                            full.names = FALSE))
  mclapply(files, create_count_overlaps_file,
           cell_number_filter=cell_number_filter,
           metadata=metadata,
           interval_ranges=interval.ranges,
           chain=ch,
           mc.cores=8)
} else if (dataset == "shendure") {
  metadata_Shendure = read.table("raw_dir/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                                  sep="\t",
                                  header=TRUE)
  # colnames(metadata_Shendure)[1] = "cell_barcode"
  # colnames(metadata_Shendure)[2] = "sample"
  # colnames(metadata_Shendure)[grep("tss", colnames(metadata_Shendure))] = "tsse"
  files_Shendure = setdiff(list.files("raw_dir/bed_files/JShendure_scATAC/"), 
                           list.dirs("raw_dir/bed_files/JShendure_scATAC/", 
                                     recursive = FALSE, 
                                     full.names = FALSE))
  lapply(files_Shendure, create_count_overlaps_file_shendure,
         cell_number_filter=cell_number_filter,
         metadata=metadata_Shendure,
         interval_ranges=interval.ranges,
         chain=ch)
} else if (dataset == "tsankov") {
  metadata_tsankov_proximal = 
    read.csv("raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv")
  metadata_tsankov_distal = 
    read.csv("raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv")
  
  files_Tsankov_distal = list.files("raw_dir/bed_files/Tsankov_scATAC/", 
                                    pattern=".*distal.*")
  files_Tsankov_proximal = list.files("raw_dir/bed_files/Tsankov_scATAC/", 
                                      pattern=".*proximal*")
  mclapply(files_Tsankov_proximal, 
           create_count_overlaps_file_tsankov,
           cell_number_filter=cell_number_filter,
           metadata=metadata_tsankov_proximal,
           interval_ranges=interval.ranges,
           mc.cores = 4)
  
  mclapply(files_Tsankov_distal,
            create_count_overlaps_file_tsankov,
            cell_number_filter=cell_number_filter,
            metadata=metadata_tsankov_distal,
            interval_ranges=interval.ranges,
            mc.cores=4)
}


#                mc.cores = 1)
# lapply(files, create_count_overlaps_file,
#        cell_number_filter=CELL_NUMBER_FILTER,
#        metadata=metadata,
#        interval_ranges=interval.ranges,
#        chain=ch)

# lapply(files_Shendure, create_count_overlaps_file_shendure,
#         cell_number_filter=CELL_NUMBER_FILTER,
#         metadata=metadata_Shendure,
#         interval_ranges=interval.ranges)
#         mc.cores = 1)

# mclapply(files_Tsankov_proximal, create_count_overlaps_file_tsankov,
#         cell_number_filter=CELL_NUMBER_FILTER,
#         metadata=metadata_tsankov_proximal,
#         interval_ranges=interval.ranges,
#         mc.cores = 4)

# lapply(files_Tsankov_proximal, create_count_overlaps_file_tsankov,
#        cell_number_filter=CELL_NUMBER_FILTER,
#        metadata=metadata_tsankov_proximal,
#        interval_ranges=interval.ranges,
#        chain=ch,
#        top_tsse_fragment_count_subsampling_range = c(1000, 3853, 14849, 57225,
#                                                      220520))

# mclapply(files_Tsankov_distal,
#           create_count_overlaps_file_tsankov,
#           cell_number_filter=CELL_NUMBER_FILTER,
#           metadata=metadata_tsankov_distal,
#           interval_ranges=interval.ranges,
#           mc.cores=4)

# top_tsse_fragment_count_range = c(1000, 3853, 14849, 57225,
#                                              220520, 849788, 3274710,
#                                              12619294)
# 
# bing_ren_lung_files = setdiff(list.files("raw_dir/bed_files/",
#                                          pattern=".*lung.*"),
#                               list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                         full.names = FALSE))
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_lung_files,
#                                               CELL_NUMBER_FILTER,
#                                               metadata,
#                                               interval.ranges,
#                                               ch,
#                                               top_tsse_fragment_count_range,
#                                               c("Alveolar Type 2 (AT2) Cell"))

# shendure_lung_files  = setdiff(list.files("raw_dir/bed_files/JShendure_scATAC/",
#                                           pattern=".*lung.*"),
#                                list.dirs("raw_dir/bed_files/JShendure_scATAC",
#                                          recursive = FALSE,
#                                          full.names = FALSE))
# create_tsse_filtered_count_overlaps_per_tissue(shendure_lung_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata_Shendure,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                c("Lung Bronchiolar and alveolar epithelial cells",
#                                                  ))

# create_tsse_filtered_count_overlaps_per_tissue(files_Tsankov_distal,
#                                                CELL_NUMBER_FILTER,
#                                                metadata_tsankov_distal,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                c("Distal Lung AT2"))



# #### Melanoma ####
# bing_ren_skin_files = setdiff(list.files("raw_dir/bed_files/",
#                                          pattern=".*skin_SM*"),
#                               list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                         full.names = FALSE))
# 
# top_tsse_fragment_count_range =  c(1000, 10000, 50000, 100000, 150000, 250000,
#                                    300000, 400000, 500000, 600000)
# skin_cell_types = c("Melanocyte",
#                     "Fibroblast (Epithelial)")
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_skin_files,
#                                               CELL_NUMBER_FILTER,
#                                               metadata,
#                                               interval.ranges,
#                                               ch,
#                                               top_tsse_fragment_count_range,
#                                               skin_cell_types,
#                                               "bing_ren")
# 
# 
# bing_ren_skin_sun_exposed_files = setdiff(list.files("raw_dir/bed_files/",
#                                               pattern=".*skin_sun_exposed.*"),
#                               list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                         full.names = FALSE))
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_skin_sun_exposed_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                skin_cell_types,
#                                                "bing_ren")

# #### ColoRect.AdenoCA ####
# bing_ren_colon_files = setdiff(list.files("raw_dir/bed_files/",
#                                          pattern=".*colon_transverse.*"),
#                               list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                         full.names = FALSE))
# 
# top_tsse_fragment_count_range = c(1000, 3005, 9035, 26558, 27162, 81648,
#                                   245435, 737778, 2217754, 6666549)
# 
# colon_cell_types = c("Colon Epithelial Cell 2",
#                      "Colonic Goblet Cell",
#                      "T Lymphocyte 1 (CD8+)")
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_colon_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                colon_cell_types,
#                                                "bing_ren")
# 
# bing_ren_mammary_tissue_files = setdiff(list.files("raw_dir/bed_files/",
#                                           pattern=".*mammary_tissue.*"),
#                                list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                          full.names = FALSE))
# 
# mammary_tissue_cell_types = c("Basal Epithelial (Mammary)",
#                               "Mammary Luminal Epithelial Cell 1")
# 
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_mammary_tissue_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                mammary_tissue_cell_types,
#                                                "bing_ren")

# #### Breast.AdenoCA ####
# bing_ren_mammary_tissue_files = setdiff(list.files("raw_dir/bed_files/",
#                                           pattern=".*mammary_tissue.*"),
#                                list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                          full.names = FALSE))
# 
# top_tsse_fragment_count_range = c(1000, 3005, 9035, 26558, 27162, 81648,
#                                   245435, 737778, 2217754, 6666549, 20039584)
# 
# mammary_tissue_cell_types = c("Basal Epithelial (Mammary)",
#                               "Mammary Luminal Epithelial Cell 1",
#                               "Mammary Luminal Epithelial Cell 2")
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_mammary_tissue_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                mammary_tissue_cell_types,
#                                                "bing_ren")
# 
# bing_ren_heart_atrial_appendage_tissue_files = setdiff(list.files("raw_dir/bed_files/",
#                                                    pattern=".*heart_atrial_appendage.*"),
#                                         list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                                   full.names = FALSE))
# 
# heart_atrial_appendage_cell_types = c("Cardiac Pericyte 1")
# 
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_heart_atrial_appendage_tissue_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                heart_atrial_appendage_cell_types,
#                                                "bing_ren")
# bing_ren_esophagus_mucosa_tissue_files = setdiff(list.files("raw_dir/bed_files/",
#                                                                   pattern=".*esophagus_mucosa.*"),
#                                                        list.dirs("raw_dir/bed_files", recursive = FALSE,
#                                                                  full.names = FALSE))
# 
# esophagus_mucosa_cell_types = c("Airway Goblet Cell")
# 
# 
# create_tsse_filtered_count_overlaps_per_tissue(bing_ren_esophagus_mucosa_tissue_files,
#                                                CELL_NUMBER_FILTER,
#                                                metadata,
#                                                interval.ranges,
#                                                ch,
#                                                top_tsse_fragment_count_range,
#                                                esophagus_mucosa_cell_types,
#                                                "bing_ren")
