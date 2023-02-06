library(parallel)
library(tidyverse)
library(rtracklayer)
library(optparse)
library(rprojroot)

# setwd((find_rstudio_root_file()))
source("/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/scripts/count_overlaps_utils.R")

option_list <- list( 
  make_option("--cell_types", type="character"),
  make_option("--dataset", type="character"),
  make_option("--top_tsse_fragment_count_range", type="character"),
  make_option("--files_pattern", type="character"),
  make_option("--cores", type="integer")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,600000",
#          "--dataset=bing_ren",
#          "--cell_types=*",
#           "--files_pattern=lung",
#          "--cores=2"))

# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--top_tsse_fragment_count_range=1000",
#          "--dataset=bing_ren",
#          "--cell_types=Colon Epithelial Cell 2,Colon Epithelial Cell 1",
#           "--files_pattern=colon_transverse",
#          "--cores=1"))

#args = parse_args(OptionParser(option_list=option_list), args=
#       c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#         "--dataset=tsankov",
#         "--cell_types=AT2",
#          "--files_pattern=RPL",
#         "--cores=1"))

# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#          "--dataset=bing_ren",
#          "--cell_types=Airway Goblet Cell,Esophageal Epithelial Cell",
#           "--files_pattern=esophagus_mucosa",
#          "--cores=1"))

# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#          "--dataset=bing_ren",
#          "--cell_types=Airway Goblet Cell,Esophageal Epithelial Cell",
#           "--files_pattern=esophagus_mucosa",
#          "--cores=1"))

# args = parse_args(OptionParser(option_list=option_list), args=
#        c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,600000",
#          "--dataset=shendure",
#          "--cell_types=Astrocytes/Oligodendrocytes,Astrocytes",
#           "--files_pattern=brain",
#          "--cores=2"))
# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,600000,1000000,5000000,10000000,14388810",
#                      "--dataset=shendure",
#                      "--cell_types=Goblet cells",
#                      "--files_pattern=stomach",
#                      "--cores=1"))
# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#                      "--dataset=bing_ren",
#                      "--cell_types=Cardiac Pericyte 1,Cardiac Pericyte 2",
#                      "--files_pattern=heart_lv_SM",
#                      "--cores=1"))
# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#                      "--dataset=bing_ren",
#                      "--cell_types=Small Intestinal Enterocyte",
#                      "--files_pattern=intestine",
#                      "--cores=1"))
#args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000",
#                      "--dataset=bing_ren",
#                      "--cell_types=Cardiac Pericyte 1,Cardiac Pericyte 2,Endothelial Cell (Myocardial)",
#                      "--files_pattern=heart_lv_SM",
#                      "--cores=1"))

# args = parse_args(OptionParser(option_list=option_list), args=
#                    c("--top_tsse_fragment_count_range=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000",
#                      "--dataset=tsankov",
#                      "--cell_types=AT1;AT2;B_cells;Ciliated;Endothelial;Fibroblasts;Immune;Mesothelium;Secretory;SmoothMuscle",
#                      "--files_pattern=RPL",
#                      "--cores=1"))

args = parse_args(OptionParser(option_list=option_list), args=
                   c("--top_tsse_fragment_count_range=100000,150000,200000,300000,400000,500000,1000000",
                     "--dataset=Tsankov",
                     "--cell_types=Basal;Ciliated;Secretory;Myeloid;Endothelial;Neuroendocrine;B.cells;Ionocytes;T.NK.cells;Stromal;Tuft.like;Sec-Ciliated",
                     "--files_pattern=RPL",
                     "--cores=1"))
top_tsse_fragment_count_range = as.integer(unlist(strsplit(
                                           args$top_tsse_fragment_count_range, 
                                           split = ",")))
dataset = args$dataset
cell_types = unlist(strsplit(args$cell_types, split = ";"))
files_pattern = args$files_pattern
cores = args$cores

get_fragments_from_top_cells <- function(i, 
                                         count_filtered_metadata_with_fragment_counts_per_cell_type, 
                                         fragments, sample_names) {
  top_barcodes = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type,
                        function(x) x[["cell_barcode"]])
  samples = lapply(count_filtered_metadata_with_fragment_counts_per_cell_type, 
                   function(x) x[["sample"]])
  samples_per_cell_type = samples[[i]]
  
  sample_idx = unlist(lapply(samples_per_cell_type, match, sample_names))
  fragments_from_relevant_barcodes = lapply(seq_along(sample_idx), 
                                            get_fragments_by_cell_barcode, 
                                            sample_idx, 
                                            fragments, 
                                            top_barcodes[[i]])
  fragments_from_relevant_barcodes = do.call(c, 
                                             fragments_from_relevant_barcodes)
  return(fragments_from_relevant_barcodes)
}

create_tsse_filtered_count_overlaps_helper <- function(metadata_with_fragment_counts,
                                                       cell_types_to_consider,
                                                       migrated_fragments,
                                                       sample_names,
                                                       interval_ranges,
                                                       count, 
                                                       cores) {
  metadata_with_fragment_counts = metadata_with_fragment_counts %>% 
                                  filter(cell_type %in% cell_types_to_consider) %>%
                                  group_by(cell_type) %>%
                                  arrange(desc(tsse), .by_group = TRUE) %>%
                                  mutate(frag_counts_cumsum = 
                                           cumsum(frag_counts))
  metadata_with_fragment_counts_per_cell_type = 
    group_split(metadata_with_fragment_counts)
  cell_types = unlist(lapply(metadata_with_fragment_counts_per_cell_type, 
                             function(x) x[["cell_type"]][1]))
  debug_message = paste0("Creating count overlaps for num fragments = ",
                         count)
  print(debug_message)
  count_filtered_metadata_with_fragment_counts_per_cell_type =
    mclapply(metadata_with_fragment_counts_per_cell_type,
             filter_metadata_by_fragment_count, count,
             mc.cores=cores)
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
                                      sample_names, mc.cores=cores)
  set.seed(42)
  fragments_from_top_cells = mclapply(fragments_from_top_cells, sample, count,
                                      mc.cores=cores)
  count_overlaps = mclapply(fragments_from_top_cells,
                            function(x) countOverlaps(interval_ranges, x),
                            mc.cores=cores)
  count_overlaps = as_tibble(count_overlaps, 
                             .name_repair="universal")
  colnames(count_overlaps) = cell_types
  return(count_overlaps)
}

create_tsse_filtered_count_overlaps_per_tissue <- function(files,
                                                           metadata,
                                                           interval_ranges, 
                                                           chain, 
                                                           top_tsse_fragment_count_range=NA, 
                                                           cell_types_to_consider=NA,
                                                           dataset="bing_ren",
                                                           cores) {
  
  if (dataset == "Bingren") {
    filepaths = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                      "data", "bed_files", "bingren_scATAC", files, sep="/")
  }
  else if (dataset == "Shendure") {
    filepaths = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                      "data", "bed_files", "JShendure_scATAC", 
                      files, sep="/")
  }
  else {
    filepaths = paste("/broad", "hptmp", "bgiotti", "BingRen_scATAC_atlas", 
                      "data", "bed_files", "Tsankov_scATAC", 
                      files, sep="/")
    print(filepaths)
  }

  print("Importing BED files...")
  fragments = mclapply(filepaths, import, format="bed", mc.cores=cores)

  if (dataset == "Bingren") {
    print("Migrating BED files...")
    migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
                                  mc.cores=cores)
    sample_names = unlist(lapply(files, get_sample_name))
  }
  else if (dataset == "Shendure") {
    migrated_fragments = fragments
    sample_names = unlist(lapply(files, get_sample_name_shendure))
  }
  else {
    print("Migrating BED files...")
    migrated_fragments = mclapply(fragments, migrate_bed_file_to_hg37, ch, 
                                  mc.cores=cores)
    sample_names = unlist(lapply(files, get_sample_name_tsankov))
  }
  # if (dataset == "bing_ren" || dataset == "shendure") {
  filtered_metadatas = mclapply(sample_names, filter_metadata_by_sample_name, 
                                metadata, mc.cores=cores)
  # }
  # else if (dataset == "tsankov") {
  #   # No need to filter anything out because metadata for Tsankov is per
  #   # tissue/sample type
  #   if (length(files) == 4) {
  #     filtered_metadatas = list(metadata, metadata, metadata, metadata)
  #   }
  #   else if (length(files) == 3) {
  #     filtered_metadatas = list(metadata, metadata, metadata)
  #   }
  # }
  
  if (dataset == "Bingren") {
    sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
                                            get_sample_barcodes_in_metadata,
                                            "Bingren")
    filtered_metadatas = lapply(seq_along(filtered_metadatas),
                                add_cell_barcodes_to_metadata,
                                filtered_metadatas,
                                sample_barcodes_in_metadatas)
  }
  else if (dataset == "Shendure") {
    sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
                                            function(x) x[["cell_barcode"]])
  }
  else {
    sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
                                            get_sample_barcodes_in_metadata,
                                            "Tsankov")
    sample_barcodes_in_metadatas = lapply(sample_barcodes_in_metadatas, 
                                          substr, 1, 16)
    fix_barcodes <- function(i, metadatas, sample_barcodes_in_metadatas) {
      metadatas[[i]]["cell_barcode"] = sample_barcodes_in_metadatas[[i]]
      return(metadatas[[i]])
    }
    
    filtered_metadatas = mclapply(seq_along(filtered_metadatas), fix_barcodes, 
                                  filtered_metadatas, sample_barcodes_in_metadatas,
                                  mc.cores=cores)
    remove_name_suffix <- function(granges) {
      granges$name = substr(granges$name, 1, 16)
      return(granges)
    }
    migrated_fragments = mclapply(migrated_fragments, remove_name_suffix, 
                                  mc.cores=cores)
    # sample$name = substr(sample$name, 1, 16)
    
    # sample_barcodes_in_metadatas = mclapply(filtered_metadatas, 
    #                                         function(x) x[["cell_barcode"]])
  }
  migrated_fragments <- mclapply(seq_along(migrated_fragments),
                                 filter_samples_to_contain_only_cells_in_metadata,
                                 migrated_fragments,
                                 sample_barcodes_in_metadatas,
                                 mc.cores=cores)
  migrated_fragments <- mclapply(migrated_fragments,
                                 function(x) x[!is.na(findOverlaps(x, interval_ranges,
                                                          select = "arbitrary"))],
                                 mc.cores=cores)
  print("Counting fragments per cell")
  fragment_counts_per_sample = mclapply(migrated_fragments, count_fragments_per_cell,
                                        mc.cores=cores)
  
  idx = 1
  metadata_with_fragment_counts = tibble()
  for (metadata in filtered_metadatas) {
    metadata["frag_counts"] = as.integer(fragment_counts_per_sample[[idx]][match(sample_barcodes_in_metadatas[[idx]], 
                                                                names(fragment_counts_per_sample[[idx]]))])
    metadata_with_fragment_counts = bind_rows(metadata_with_fragment_counts, 
                                              as_tibble(metadata))
    # filtered_metadatas[idx] = metadata
    idx = idx + 1
  }
  
  # PARALLELIZED
  # mclapply(top_tsse_fragment_count_range, 
  #          create_count_overlaps_file_per_tsse_count, 
  #          metadata_with_fragment_counts_per_cell_type,
  #          migrated_fragments,
  #          sample_names, interval_ranges,
  #          cell_types_to_consider,
  #          mc.cores=8)
  
  metadata_with_fragment_counts = metadata_with_fragment_counts %>% 
                                  filter(cell_type %in% cell_types_to_consider)
  
  if (dataset == "Bingren") {
    tissue_name = get_tissue_name(files[1])
  }
  else if (dataset == "Shendure") {
    tissue_name = get_tissue_name_shendure(files[1])
  }
  else {
    if (grepl("RPL", files[1])) {
      print("distal_lung")
      tissue_name = "distal_lung"
    }
    else {
      print("proximal lung")
      tissue_name = "proximal_lung"
    }
  }
  
  for (count in top_tsse_fragment_count_range) {
    filename = paste("count_overlaps", "frag_count_filter", count, sep="_")
    filename = paste(filename, "rds", sep=".")
    filepath = paste("/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/count_overlap_data/tsse_filtered",
                     dataset, tissue_name, sep="/")
    dir.create(filepath, recursive = T)
    filepath = paste(filepath, filename, sep="/")
    if (file.exists(filepath)) {
      count_overlaps = readRDS(filepath)
      already_existing_cell_types = colnames(count_overlaps)
      cell_types_to_consider_temp = cell_types_to_consider[!(cell_types_to_consider %in% 
                                                               already_existing_cell_types)]
      if (length(cell_types_to_consider_temp != 0)) {
        count_overlaps_from_missing_cells = 
          create_tsse_filtered_count_overlaps_helper(metadata_with_fragment_counts,
                                                     cell_types_to_consider_temp,
                                                     migrated_fragments,
                                                     sample_names, 
                                                     interval_ranges,
                                                     count,
                                                     cores)
        # for (co in count_overlaps_from_missing_cells) {
        if (length(count_overlaps_from_missing_cells) > 0) {
          if (length(count_overlaps) > 0) {
            count_overlaps = cbind(count_overlaps, 
                                   count_overlaps_from_missing_cells)
          }
          else {
            count_overlaps = count_overlaps_from_missing_cells
          }
          saveRDS(count_overlaps, filepath)
        }
        # }
      }
    }
    else {
      count_overlaps = create_tsse_filtered_count_overlaps_helper(metadata_with_fragment_counts,
                                                                  cell_types_to_consider,
                                                                  migrated_fragments,
                                                                  sample_names, 
                                                                  interval_ranges,
                                                                  count,
                                                                  cores)
      saveRDS(count_overlaps, filepath)
    }
  }
}

load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
# hg38_path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/hg38ToHg19.over.chain")

if (dataset == "Bingren") {
  metadata = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/GSE184462_metadata.tsv", sep="\t", 
                        header=TRUE)
  colnames(metadata)[grep("cell.type", colnames(metadata))] = "cell_type"
  files = setdiff(list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/", pattern=files_pattern,
                             full_names=TRUE),
                  list.dirs("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files", recursive = FALSE, full.names = 
                             TRUE))
  create_tsse_filtered_count_overlaps_per_tissue(files,
                                                 metadata,
                                                 interval.ranges,
                                                 ch,
                                                 top_tsse_fragment_count_range,
                                                 cell_types,
                                                 dataset,
                                                 cores)
} else if (dataset == "Shendure") {
  metadata = read.table("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                        sep="\t",
                        header=TRUE)
  colnames(metadata)[1] = "cell_barcode"
  colnames(metadata)[2] = "sample"
  colnames(metadata)[grep("frit", colnames(metadata))] = "tsse"
  
  files = setdiff(list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/JShendure_scATAC/", 
                             pattern=files_pattern, full.names = TRUE),
                  list.dirs("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/JShendure_scATAC/", 
                            recursive = FALSE, full.names = TRUE))
  create_tsse_filtered_count_overlaps_per_tissue(files,
                                                 metadata,
                                                 interval.ranges,
                                                 ch,
                                                 top_tsse_fragment_count_range,
                                                 cell_types,
                                                 dataset,
                                                 cores)
} else if (dataset == "Tsankov") {
  metadata_tsankov_proximal = 
    read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/tsankov_lung_proximal_barcode_annotation.csv")
  metadata_tsankov_distal = 
    read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/tsankov_lung_distal_barcode_annotation.csv")
  
  colnames(metadata_tsankov_proximal)[1] = "cell_barcode"
  colnames(metadata_tsankov_proximal)[2] = "sample"
  colnames(metadata_tsankov_proximal)[grep("TSSEnrichment", 
                                 colnames(metadata_tsankov_proximal))] = "tsse"
  colnames(metadata_tsankov_proximal)[grep("celltypes", 
                            colnames(metadata_tsankov_proximal))] = "cell_type"
  
  colnames(metadata_tsankov_distal)[1] = "cell_barcode"
  colnames(metadata_tsankov_distal)[2] = "sample"
  colnames(metadata_tsankov_distal)[grep("TSSEnrichment", 
                                  colnames(metadata_tsankov_distal))] = "tsse"
  colnames(metadata_tsankov_distal)[grep("celltypes", 
                              colnames(metadata_tsankov_distal))] = "cell_type"
  
  files_Tsankov_proximal = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/Tsankov_scATAC/", 
                                      pattern="IC.*[.]gz")
  files_Tsankov_distal = list.files("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/bed_files/Tsankov_scATAC/", 
                                    pattern="RPL.*[.]gz")
  if (files_pattern == "RPL") {
    print("Tsankov RPL")
    files_Tsankov_distal = c("RPL_280_neg_1_fragments.tsv.gz", "RPL_280_neg_2_fragments.tsv.gz")
    create_tsse_filtered_count_overlaps_per_tissue(files=files_Tsankov_distal,
                                                   metadata=metadata_tsankov_distal,
                                                   interval_ranges=interval.ranges,
                                                   chain=ch,
                                                   top_tsse_fragment_count_range,
                                                   cell_types,
                                                   dataset,
                                                   cores)
      
  }
  else {
    print("Tsankov IC")
    create_tsse_filtered_count_overlaps_per_tissue(files=files_Tsankov_proximal,
                                                   metadata=metadata_tsankov_proximal,
                                                   interval_ranges=interval.ranges,
                                                   chain=ch,
                                                   top_tsse_fragment_count_range,
                                                   cell_types_to_consider=cell_types,
                                                   dataset,
                                                   cores)
  }
}
