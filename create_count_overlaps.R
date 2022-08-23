library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(exomeCopy)
library(parallel)
library(liftOver)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied", call.=FALSE)
} else {
  CELL_NUMBER_FILTER = as.integer(args[1])
} 

# add_cumulative_num_fragments_per_cell_to_metadata <- function(filtered_metadata,
#                                                               sample) {
#   sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(
#                                             filtered_metadata, "cellID", "\\+")
#   filtered_metadata["cell_barcode"] = sample_barcodes_in_metadata
#   frag_counts_per_cell = table(sample$name)
#   match(filtered_metadata[["cell_barcode"]], names(frag_counts_per_cell))
# }

# match_barcodes_to_fragment_counts <- function(i, barcodes, fragment_counts) {
#   fragment_counts[match(barcodes[i], names(fragment_counts))]
# }

add_cell_barcodes_to_metadata <- function(i, metadatas, barcodes) {
  metadatas[[i]]["cell_barcode"] = barcodes[i]
  return(metadatas[[i]])
}

add_fragment_counts_per_cell_to_metadata <- function(metadata, barcodes,
                                                     fragment_counts_per_sample) {
  metadata["counts"] = fragment_counts_per_sample[match(barcodes, 
                                             names(fragment_counts_per_sample))]
  return(metadata)
  # metadatas[[i]]["counts"] = 
  #   fragment_counts_per_sample[[i]][match(barcodes[[i]], 
  #                                         names(fragment_counts_per_sample[[i]]))]
}

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

count_fragments_per_cell <- function(fragments) {
  return(table(fragments$name))
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
    filtered_metadata = metadata[metadata$sample %like% sample_name, ]
    sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(filtered_metadata,
                                                                  "cellID",
                                                                  "\\+")
    sample <- filter_sample_to_contain_only_cells_in_metadata(sample,
                                                    sample_barcodes_in_metadata)
    sample <- get_sample_cell_types(sample, filtered_metadata)
    
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
    sample <- get_sample_cell_types_shendure(sample, sample_name, metadata)
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

create_tsse_filtered_count_overlaps_helper <- function(metadata_with_fragment_counts,
                                                       cell_types_to_consider,
                                                       migrated_fragments,
                                                       sample_names,
                                                       interval_ranges,
                                                       count) {
  metadata_with_fragment_counts = metadata_with_fragment_counts %>% 
                                  filter(cell_type %in% cell_types_to_consider) %>%
                                  group_by(cell_type) %>%
                                  arrange(desc(tsse), .by_group = TRUE) %>%
                                  mutate(frag_counts_cumsum = 
                                           cumsum(frag_counts))
  metadata_with_fragment_counts_per_cell_type = group_split(metadata_with_fragment_counts)
  cell_types = unlist(lapply(metadata_with_fragment_counts_per_cell_type, 
                      function(x) x[["cell_type"]][1]))
  debug_message = paste0("Creating count overlaps for num fragments = ",
                         count)
  print(debug_message)
  count_filtered_metadata_with_fragment_counts_per_cell_type =
                            mclapply(metadata_with_fragment_counts_per_cell_type,
                                     filter_metadata_by_fragment_count, count,
                                     mc.cores=8)
  cell_types_filter = !unlist(lapply(count_filtered_metadata_with_fragment_counts_per_cell_type,
                                    is.null))
  cell_types = cell_types[cell_types_filter]
  # cell_types_to_consider = cell_types_to_consider[cell_types_filter]
  count_filtered_metadata_with_fragment_counts_per_cell_type = 
    count_filtered_metadata_with_fragment_counts_per_cell_type[cell_types_filter]
  # top_barcodes = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type,
  #                       function(x) x[["cell_barcode"]])
  fragments_from_top_cells = mclapply(seq_along(count_filtered_metadata_with_fragment_counts_per_cell_type),
                                      get_fragments_from_top_cells,
                                      count_filtered_metadata_with_fragment_counts_per_cell_type,
                                      migrated_fragments,
                                      sample_names, mc.cores=8)
                              
  fragments_from_top_cells = mclapply(fragments_from_top_cells, sample, count,
                                      mc.cores=8)
  count_overlaps = mclapply(fragments_from_top_cells,
                            function(x) countOverlaps(interval_ranges, x),
                            mc.cores=8)
  count_overlaps = as_tibble(count_overlaps, 
                             .name_repair="universal")
  colnames(count_overlaps) = cell_types
  return(count_overlaps)
}

create_tsse_filtered_count_overlaps_per_tissue <- function(files, 
                                                           cell_number_filter, 
                                                           metadata,
                                                           interval_ranges, 
                                                           chain, 
                                                           top_tsse_fragment_count_range=NA, 
                                                           cell_types_to_consider=NA,
                                                           dataset="bing_ren") {
  
  if (dataset == "bing_ren") {
    filepaths = paste("raw_dir", "bed_files", files, sep="/")
  }
  else if (dataset == "shendure") {
    filepaths = paste("raw_dir", "bed_files", "JShendure_scATAC", 
                      files, sep="/")
  }
  else {
    filepaths = paste("raw_dir", "bed_files", "Tsankov_scATAC", 
                      files, sep="/")
  }
  # # ##### DELETE #####
  # fragments = lapply(c("raw_dir/bed_files/temp_GSM5589377_lung_SM-A8WNH_rep1_fragments.rds",
  #                      "raw_dir/bed_files/temp_GSM5589378_lung_SM-ACCPU_rep1_fragments.rds"),
  #                    readRDS)
  # # ##################
  print("Importing BED files...")
  fragments = mclapply(filepaths, import, format="bed", mc.cores=8)
  print("Migrating BED files...")
  migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
                                mc.cores=8)
  
  if (dataset == "bing_ren") {
    sample_names = unlist(lapply(files, get_sample_name))
  }
  else if (dataset == "shendure") {
    sample_names = unlist(lapply(files, get_sample_name_shendure))
  }
  
  # # ##### DELETE #####
  # sample_names = c("lung_SM-ACCPU_1", "lung_SM-A8WNH_1")
  # # ##################
  
  if (dataset == "bing_ren" || dataset == "shendure") {
    filtered_metadatas = mclapply(sample_names, filter_metadata_by_sample_name, 
                                  metadata, mc.cores=8)
  }
  else if (dataset == "tsankov") {
    filtered_metadatas = metadata
  }
  
  print("Counting fragments per cell")
  fragment_counts_per_sample = mclapply(migrated_fragments, count_fragments_per_cell,
                                        mc.cores=8)
  if (dataset == "bing_ren") {
    sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
                                            get_sample_barcodes_in_metadata,
                                            "cellID", 
                                            "\\+")
    filtered_metadatas = lapply(seq_along(filtered_metadatas),
                                add_cell_barcodes_to_metadata,
                                filtered_metadatas,
                                sample_barcodes_in_metadatas) 
  }
  else if (dataset == "shendure") {
    sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
                                            function(x) x[["cell_barcode"]])
  }


  # sum(!is.na(match(sample_barcodes_in_metadatas[[2]], 
  #                  names(fragment_counts_per_sample)[[2]])))
  
  idx = 1
  metadata_with_fragment_counts = tibble()
  for (metadata in filtered_metadatas) {
    metadata["frag_counts"] = fragment_counts_per_sample[[idx]][match(sample_barcodes_in_metadatas[[idx]], 
                                      names(fragment_counts_per_sample[[idx]]))]
    metadata_with_fragment_counts = bind_rows(metadata_with_fragment_counts, 
                                              as_tibble(metadata))
    # filtered_metadatas[idx] = metadata
    idx = idx + 1
  }
  # sample_barcodes_in_metadatas = get_sample_barcodes_in_metadata(
  #                                       filtered_metadata,
  #                                       "cellID", "\\+")
  
  # metadata_with_fragment_counts = metadata_with_fragment_counts %>%
  #                                 filter(cell_type %in% cell_types_to_consider) %>%
  #                                 group_by(cell_type) %>%
  #                                 arrange(desc(tsse), .by_group = TRUE) %>%
  #                                 mutate(frag_counts_cumsum =
  #                                          cumsum(frag_counts))
  # metadata_with_fragment_counts_per_cell_type = group_split(metadata_with_fragment_counts)
  
  # metadata_with_fragment_counts_per_cell_type["frag_counts_cumsum"] = 
  #   cumsum(metadata_with_fragment_counts_per_cell_type$counts)
  
  # PARALLELIZED
  # mclapply(top_tsse_fragment_count_range, 
  #          create_count_overlaps_file_per_tsse_count, 
  #          metadata_with_fragment_counts_per_cell_type,
  #          migrated_fragments,
  #          sample_names, interval_ranges,
  #          cell_types_to_consider, 
  #          cell_number_filter,
  #          mc.cores=8)
  metadata_with_fragment_counts = metadata_with_fragment_counts %>% 
                                  filter(cell_type %in% cell_types_to_consider)
  
  for (count in top_tsse_fragment_count_range) {
    if (dataset == "bing_ren") {
      filename = paste("count_filter", cell_number_filter,
                     "count_overlaps", "frag_count_filter", count, sep="_")
      tissue_name = get_tissue_name(files[1])
    }
    else if (dataset == "shendure") {
      filename = paste("shendure_count_filter", cell_number_filter,
                       "count_overlaps", "frag_count_filter", count, sep="_")
      tissue_name = get_tissue_name_shendure(files[1])
    }
    filename = paste(filename, "rds", sep=".")
    filepath = paste("processed_data/count_overlap_data/tsse_filtered",
                     tissue_name, sep="/")
    dir.create(filepath)
    filepath = paste(filepath, filename, sep="/")
    if (file.exists(filepath)) {
      count_overlaps = readRDS(filepath)
      already_existing_cell_types = colnames(count_overlaps)
      cell_types_to_consider = cell_types_to_consider[!(cell_types_to_consider %in% 
                                                      already_existing_cell_types)]
      if (length(cell_types_to_consider != 0)) {
        count_overlaps_from_missing_cells = 
          create_tsse_filtered_count_overlaps_helper(metadata_with_fragment_counts,
                                                     cell_types_to_consider,
                                                     migrated_fragments,
                                                     sample_names, 
                                                     interval_ranges,
                                                     count)
        count_overlaps = cbind(count_overlaps, 
                               count_overlaps_from_missing_cells)
        saveRDS(count_overlaps, filepath)
      }
    }
    else {
      count_overlaps = create_tsse_filtered_count_overlaps_helper(metadata_with_fragment_counts,
                                                                  cell_types_to_consider,
                                                                  migrated_fragments,
                                                                  sample_names, 
                                                                  interval_ranges,
                                                                  count)
      saveRDS(count_overlaps, filepath)
    }
  }
}

  # # 
  # # l = match(sample_barcodes_in_metadatas, fragments$name)
  # filtered_fragments = c()
  # idx = 1
  # for (sample in fragments) {
  #   filtered_fragments = append(filtered_fragments,
  #                               filter_sample_to_contain_only_cells_in_metadata(sample,
  #                               sample_barcodes_in_metadatas[[idx]]))
  #   idx = idx + 1
  # }
  # 
  # a <- mclapply(seq_along(fragments), 
  #                       filter_sample_to_contain_only_cells_in_metadata,
  #                       fragments,
  #                       sample_barcodes_in_metadatas)  
  # 
  # fragments = lapply(seq_along(fragments), get_sample_cell_types, 
  #                    fragments, sample_barcodes_in_metadatas,
  #                    filtered_metadatas)
  # 
  # 
  # fragments = as_tibble(fragments) %>%
  #             filter(cell.type %in% cell_types_to_consider) %>%
  #             groupby(cell.type) %>%
  #             arrange(desc())



filter_metadata_by_fragment_count <- function(metadata, min_frag_counts) {
  if ((tail(metadata["frag_counts_cumsum"], n=1) %>% pull) < min_frag_counts) {
      return(NULL)
  }
  metadata_cut_index_minus_one = sum(metadata["frag_counts_cumsum"] < 
                                       min_frag_counts)
  metadata_cut_index = metadata_cut_index_minus_one + 1
  metadata <- metadata[1:metadata_cut_index, ]
  return(metadata)
}

filter_metadata_by_sample_name <- function(sample_name, metadata) {
  return(metadata[metadata$sample %like% sample_name, ])
}

filter_sample_by_cell_number <- function(sample, 
                                         counts_per_cell_type, 
                                         cell_number_filter) {
  cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= 
                                          cell_number_filter, ]$cell_type
  sample = sample %>%                                             
           filter(cell_type %in% cell_type_keep)
  return(sample)
}

# filter_sample_by_tss <- function(sample, subsampling_frag_count, metadata) {
#   relevant_cell_types = unique(sample[["cell_type"]])
#   filtered_metadata = metadata %>% 
#                       filter(cell.type %in% relevant_cell_types)
#   filtered_metadata = filtered_metadata %>% 
#                       group_by(cell.type) %>%
#                       arrange(desc(tsse))
#   filtered_metadata = 
#             add_cumulative_num_fragments_per_cell_to_metadata(filtered_metadata,
#                                                               sample)
# }

filter_sample_to_contain_only_cells_in_metadata <- function(sample,
                                                            sample_barcodes_in_metadata) {
  sample = sample[sample$name %in% sample_barcodes_in_metadata]
  return(sample)
}

get_fragments_by_cell_barcode <- function(i, sample_idx, fragments, 
                                          top_barcodes) {
  return(fragments[[sample_idx[i]]][fragments[[sample_idx[i]]]$name %in% 
                                    top_barcodes[i]])
}

get_fragments_from_top_cells <- function(i, 
                                         count_filtered_metadata_with_fragment_counts_per_cell_type, 
                                         fragments, sample_names) {
  top_barcodes = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type,
                        function(x) x[["cell_barcode"]])
  samples = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type, 
                   function(x) x[["sample"]])
  # count = 1
  # frags_per_cell_type = c()
  samples_per_cell_type = samples[[i]]
  
  sample_idx = unlist(lapply(samples_per_cell_type, match, sample_names))
  fragments_from_relevant_barcodes = lapply(seq_along(sample_idx), 
                                            get_fragments_by_cell_barcode, 
                                            sample_idx, 
                                            fragments, 
                                            top_barcodes[[i]])
  fragments_from_relevant_barcodes = do.call(c, 
                                             fragments_from_relevant_barcodes)
  # frags_per_cell_type[[count]] = fragments_from_relevant_barcodes
  
  # for (samples_per_cell_type in samples) {
  #   sample_idx = unlist(lapply(samples_per_cell_type, match, sample_names))
  #   fragments_from_relevant_barcodes = lapply(seq_along(sample_idx), 
  #                                             get_fragments_by_cell_barcode, 
  #                                             sample_idx, 
  #                                             fragments, 
  #                                             top_barcodes[[count]])
  #   fragments_from_relevant_barcodes = do.call(c, 
  #                                              fragments_from_relevant_barcodes)
  #   frags_per_cell_type[[count]] = fragments_from_relevant_barcodes
  #   
  #   count = count + 1
  # }
  return(fragments_from_relevant_barcodes)
}

get_num_cells_per_sample <- function(sample) {
  counts_per_cell_type = sample %>%                               
                         group_by(cell_type) %>%                 
                         summarize(n_cells = n_distinct(name))
  return(counts_per_cell_type)
}

get_sample_barcodes_in_metadata <- function(metadata, cellID_col_name,
                                            separator_char) {
  sample_barcodes_in_metadata = unlist(lapply(metadata[[cellID_col_name]],
                                              function(x) unlist(strsplit(x,
                                              separator_char))[2]))
  return(sample_barcodes_in_metadata)
}

get_sample_cell_types <- function(i, fragments, sample_barcodes_in_metadatas,
                                  filtered_metadatas) {
  sample_idx_in_metadata = match(fragments[[i]]$name, 
                                 sample_barcodes_in_metadatas[[i]])
  fragments[[i]]$cell_type = filtered_metadatas[[i]][unlist(sample_idx_in_metadata), 
                                           "cell.type"]
  return(fragments)
}

get_sample_cell_types_shendure <- function(sample, sample_name, metadata) {
  filtered_metadata = metadata[metadata$sample %like% sample_name, ]
  sample_barcodes_in_metadata = filtered_metadata[["cell"]]                  
  sample = sample[sample$name %in% sample_barcodes_in_metadata]   
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                       "cell_type"]
  return(as.tibble(sample))
}

get_sample_cell_types_tsankov <- function(sample, metadata) {
  # sample_barcodes_in_metadata = unlist(lapply(metadata[["X"]],
  #                                             function(x) unlist(strsplit(x,
  #                                                                   '#'))[2]))
  sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(metadata, 
                                                                "X",
                                                                "#")
  sample = sample[sample$name %in% sample_barcodes_in_metadata]
  sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
  sample$cell_type = metadata[unlist(sample_idx_in_metadata),
                              "celltypes"]
  return(as.tibble(sample))
}

get_sample_name <- function(file) {
  # if (grepl("frontal_cortex", file)) {
  #   sample_name = ""
  # } 
  # else {
  #   sample_name = str_extract(file, pattern="_.*_SM")              
  #   sample_name = substr(sample_name, 2, nchar(sample_name) - 3)    
  # }
  sample_name = str_remove(file, "rep")
  sample_name = str_remove(sample_name, "_fragments.bed.gz")
  sample_name = unlist(strsplit(sample_name, split="_"))
  sample_name = paste(sample_name[2:length(sample_name)], collapse="_")
  return(sample_name)
}

get_sample_name_shendure <- function(file) {
  sample_name = unlist(strsplit(file, split="\\."))[1]
  return(sample_name)
}

get_tissue_name <- function(file) {
  if (grepl("frontal_cortex", file)) {
    tissue_name = "frontal_cortex"
  } 
  else {
    tissue_name = str_extract(file, pattern="_.*_SM")              
    tissue_name = substr(tissue_name, 2, nchar(tissue_name) - 3)    
  }
  return(tissue_name)
}

get_tissue_name_shendure <- function(file) {
  sample_name = get_sample_name_shendure(file)
  tissue_name = unlist(strsplit(sample_name, split="_"))[3]
  return(tissue_name)
}
# match_barcodes_to_fragment_counts <- function(i, barcodes, fragment_counts) {
#   fragment_counts[match(barcodes[i], names(fragment_counts))]
# }

migrate_bed_file_to_hg37 <- function(bed_sample, chain) {
  seqlevelsStyle(bed_sample) = "UCSC"  # necessary
  bed_sample = unlist (liftOver(bed_sample, ch))
  return(bed_sample)
}

hg38_path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(hg38_path)

metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", 
                      header=TRUE)
colnames(metadata)[grep("cell.type", colnames(metadata))] = "cell_type"

metadata_Shendure = read.table("raw_dir/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                               sep="\t",
                               header=TRUE)
colnames(metadata_Shendure)[2] = "sample"
colnames(metadata_Shendure)[1] = "cell_barcode"
colnames(metadata_Shendure)[grep("tss", colnames(metadata_Shendure))] = "tsse"

# metadata_tsankov_proximal = read.table("raw_dir/metadata/tsankov_lung_proximal_barcode_annotation.csv",
#                                        sep=",",
#                                        header=TRUE)
# metadata_tsankov_distal = read.table("raw_dir/metadata/tsankov_lung_distal_barcode_annotation.csv",
#                                        sep=",",
#                                        header=TRUE)
# colnames(metadata_tsankov_proximal)[2] <- "celltypes"
# colnames(metadata_tsankov_distal)[2] <- "celltypes"

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
interval.ranges

dir.create("processed_data/count_overlap_data", recursive=TRUE)                                       
dir.create("processed_data/cell_counts_per_sample")                                       

files = setdiff(list.files("raw_dir/bed_files/"), 
                list.dirs("raw_dir/bed_files", recursive = FALSE, 
                          full.names = FALSE))
files_Shendure = setdiff(list.files("raw_dir/bed_files/JShendure_scATAC/"), 
                         list.dirs("raw_dir/bed_files/JShendure_scATAC/", 
                                   recursive = FALSE, 
                                   full.names = FALSE))

files_Tsankov_distal = list.files("raw_dir/bed_files/Tsankov_scATAC/", 
                                  pattern=".*distal.*")
files_Tsankov_proximal = list.files("raw_dir/bed_files/Tsankov_scATAC/", 
                                    pattern=".*proximal*")

#mclapply(files, create_count_overlaps_file,
#                cell_number_filter=CELL_NUMBER_FILTER,
#                metadata=metadata,
#                interval_ranges=interval.ranges,
#                chain=ch, 
#		mc.cores=8)

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

top_tsse_fragment_count_range = c(1000, 3853, 14849, 57225,
                                             220520, 849788, 3274710,
                                             12619294)
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

shendure_lung_files  = setdiff(list.files("raw_dir/bed_files/JShendure_scATAC/",
                                          pattern=".*lung.*"),
                               list.dirs("raw_dir/bed_files/JShendure_scATAC",
                                         recursive = FALSE,
                                         full.names = FALSE))
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
# top_tsse_fragment_count_range = c(1000, 2983, 8901, 26558, 79237, 236409,
#                                   705340, 2104417)
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
# 
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

#### Breast.AdenoCA ####
bing_ren_mammary_tissue_files = setdiff(list.files("raw_dir/bed_files/",
                                          pattern=".*mammary_tissue.*"),
                               list.dirs("raw_dir/bed_files", recursive = FALSE,
                                         full.names = FALSE))

top_tsse_fragment_count_range = c(1000, 3005, 9035, 26558, 27162, 81648,
                                  245435, 737778, 2217754, 6666549, 20039584)

mammary_tissue_cell_types = c("Basal Epithelial (Mammary)",
                              "Mammary Luminal Epithelial Cell 1",
                              "Mammary Luminal Epithelial Cell 2")

create_tsse_filtered_count_overlaps_per_tissue(bing_ren_mammary_tissue_files,
                                               CELL_NUMBER_FILTER,
                                               metadata,
                                               interval.ranges,
                                               ch,
                                               top_tsse_fragment_count_range,
                                               mammary_tissue_cell_types,
                                               "bing_ren")

bing_ren_heart_atrial_appendage_tissue_files = setdiff(list.files("raw_dir/bed_files/",
                                                   pattern=".*heart_atrial_appendage.*"),
                                        list.dirs("raw_dir/bed_files", recursive = FALSE,
                                                  full.names = FALSE))

heart_atrial_appendage_cell_types = c("Cardiac Pericyte 1")


create_tsse_filtered_count_overlaps_per_tissue(bing_ren_heart_atrial_appendage_tissue_files,
                                               CELL_NUMBER_FILTER,
                                               metadata,
                                               interval.ranges,
                                               ch,
                                               top_tsse_fragment_count_range,
                                               heart_atrial_appendage_cell_types,
                                               "bing_ren")
bing_ren_esophagus_mucosa_tissue_files = setdiff(list.files("raw_dir/bed_files/",
                                                                  pattern=".*esophagus_mucosa.*"),
                                                       list.dirs("raw_dir/bed_files", recursive = FALSE,
                                                                 full.names = FALSE))

esophagus_mucosa_cell_types = c("Airway Goblet Cell")


create_tsse_filtered_count_overlaps_per_tissue(bing_ren_esophagus_mucosa_tissue_files,
                                               CELL_NUMBER_FILTER,
                                               metadata,
                                               interval.ranges,
                                               ch,
                                               top_tsse_fragment_count_range,
                                               esophagus_mucosa_cell_types,
                                               "bing_ren")
