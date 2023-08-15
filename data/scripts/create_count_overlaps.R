library(tools)                                                                  
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 
library(parallel)
library(optparse)
library(rtracklayer)
library(readxl)
source("count_overlaps_utils.R")

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--dataset_subsets", type="character", default=NULL),
  make_option("--cores", type="integer"),
  make_option("--annotation", type="character"),
  make_option("--which_interval_ranges", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args =
#                  c("--dataset=Rawlins_fetal_lung",
#                    "--cores=8",
#                    "--annotation=default_annotation",
#                    "--which_interval_ranges=polak"))

cores = args$cores
dataset = args$dataset
dataset_subsets = args$dataset_subsets
annotation = args$annotation
which_interval_ranges = args$which_interval_ranges

get_and_save_num_cells_per_sample <- function(sample, sample_file_name,
                                              annotation) {
  sample_file_name = paste0("cell_counts_", sample_file_name, ".rds")
  # sample_file_name = paste0("cell_counts_", file_path_sans_ext(sample_file_name, 
  #                                                              TRUE), ".rds")
  counts_per_cell_type = get_num_cells_per_sample(sample)
  path = "../processed_data/cell_counts_per_sample"
  file_path = paste(path, annotation, sample_file_name, sep="/")
  saveRDS(counts_per_cell_type, file_path)
  return(counts_per_cell_type)
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

create_count_overlaps_files <- function(file, metadata, interval_ranges, chain, 
                                        dataset, annotation, which_interval_ranges) {
  filename = get_sample_filename(file, dataset)
  if (which_interval_ranges != "polak") {
    filename = paste("interval_ranges", which_interval_ranges, filename,
                     sep = "_")
  }
  
  dirpath = paste("../processed_data/count_overlap_data", annotation, sep="/")
  filepath = paste(dirpath, filename, sep="/")
  dir.create(dirpath)
    if (!file.exists(filepath)) {
      print(paste("Processing", file, sep= " "))
      sample = import_sample(file, dataset, which_interval_ranges)
      sample_name = get_sample_name(file, dataset)
      filtered_metadata = filter_metadata_by_sample_name(sample_name, metadata)
      if (nrow(filtered_metadata) == 0) {
        print(paste(file, "has no high quality fragments", sep=" "))
        return()
      }
      
      sample_barcodes_in_metadata = get_sample_barcodes_in_metadata(filtered_metadata,
                                                                    dataset)
      
      if (dataset == "Tsankov" || dataset == "Greenleaf_brain" || 
          dataset == "Greenleaf_pbmc_bm" || dataset == "Yang_kidney" ||
          dataset == "Greenleaf_colon" || dataset == "Rawlins_fetal_lung") {
        sample$name = substr(sample$name, 1, 16)
      }
      
      # if (dataset == "Tsankov" || dataset == "Greenleaf_brain" || dataset == "Bingren") {
      #   sample = migrate_bed_file_to_hg37(sample, chain)
      # }

      # Because the function applied is used in another place over multiple 
      # samples, not just one, as is done here.
      sample <- unlist(lapply(c(1), 
                              filter_samples_to_contain_only_cells_in_metadata,
                              list(sample),
                              list(sample_barcodes_in_metadata)))
      
      sample <- get_sample_cell_types(sample[[1]], sample_barcodes_in_metadata,
                                      filtered_metadata)
      counts_per_cell_type <- get_and_save_num_cells_per_sample(sample, 
                                                             gsub(".rds", 
                                                                  filename),
                                                             annotation)
      # sample <- filter_sample_by_cell_number(sample,
      #                                        counts_per_cell_type, 
      #                                        cell_number_filter)
      if (dataset == "Yang_kidney") {
        sample$seqnames = paste0("chr", sample$seqnames)
      }
      
      count_overlaps <- compute_count_overlaps(sample, interval_ranges)
      saveRDS(count_overlaps, filepath)
    }
}
# 
# filter_sample_by_cell_number <- function(sample, 
#                                          counts_per_cell_type, 
#                                          cell_number_filter) {
#   cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= 
#                                           cell_number_filter, ]$cell_type
#   sample = sample %>%                                             
#            filter(cell_type %in% cell_type_keep)
#   return(sample)
# }

get_num_cells_per_sample <- function(sample) {
  counts_per_cell_type = sample %>%                               
                         group_by(cell_type) %>%                 
                         summarize(n_cells = n_distinct(name))
  return(counts_per_cell_type)
}

get_sample_cell_types <- function(fragments, 
                                  sample_barcodes_in_metadata,
                                  filtered_metadata) {
  sample_idx_in_metadata = match(fragments$name, 
                                 sample_barcodes_in_metadata)
  fragments$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), 
                                          "cell_type"]
  return(as.tibble(fragments))
}

if (which_interval_ranges == "polak") {
  load('../mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
} else if (which_interval_ranges == "yang") {
  load("../peak_set_yang.Rdata")
} else if (which_interval_ranges == "10kb") {
  if (file.exists("10kb_interval_ranges.Rdata")) {
    load("../10kb_interval_ranges.Rdata")
  } else {
    interval.ranges = read.table("../chr_Sorted_interval_ranges_10kb.tsv")[, 1:3]
    colnames(interval.ranges) = c("seqnames", "start","end")
    interval.ranges = GRanges(interval.ranges)
    save(interval.ranges, file="../10kb_interval_ranges.Rdata")
  }
}

dir.create("../processed_data/count_overlap_data", recursive=TRUE)   
cell_counts_dir = paste("../processed_data/cell_counts_per_sample", 
                        annotation, sep = "/")
dir.create(cell_counts_dir)
ch = import.chain("../hg38ToHg19.over.chain")

if (dataset == "Bingren") {
  if (annotation == "default_annotation") {
    metadata_bingren = read.table("../metadata/GSE184462_metadata.tsv", 
                          sep="\t",
                          header=T)
  } else if (annotation == "bingren_remove_same_celltype_indexing") {
    metadata_bingren = read.table("../metadata/bingren_remove_same_celltype_indexing.csv", 
                                  sep=",",
                                  header=T)
  }
  colnames(metadata_bingren)[grepl("cell.type", colnames(metadata_bingren))] = "cell_type"
  files_bingren = list.files("../bed_files/bingren_scATAC/migrated_to_hg19",
                              pattern=".*fragments\\.bed\\.bgz$")
  
  mclapply(files_bingren, create_count_overlaps_files,
           metadata=metadata_bingren,
           interval_ranges=interval.ranges,
           which_interval_ranges=which_interval_ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           mc.cores=cores)
} else if (dataset == "Shendure") {
  metadata_Shendure = read.table("../metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                                  sep="\t",
                                  header=TRUE)
  files_Shendure = list.files("../bed_files/JShendure_scATAC/",
                              pattern = ".*fragments\\.txt\\.gz")
  
  mclapply(files_Shendure, create_count_overlaps_files,
         metadata=metadata_Shendure,
         interval_ranges=interval.ranges,
         chain=ch,
         dataset=dataset,
         annotation=annotation,
         mc.cores=cores)
} else if (dataset == "Tsankov") {
  # metadata_tsankov_proximal = 
  #   read.csv("../metadata/tsankov_lung_proximal_barcode_annotation.csv")
  # metadata_tsankov_distal = 
  #   read.csv("../metadata/tsankov_lung_distal_barcode_annotation.csv")
  if (annotation == "Tsankov_refined") {
    metadata = read.csv("../metadata/tsankov_refined_annotation.csv")
    colnames(metadata) = c("sample", "cell_type")
    files_Tsankov_distal = list.files("../bed_files/Tsankov_scATAC/migrated_to_hg19/", 
                                      pattern="RPL.*bgz$")
    files_Tsankov_proximal = list.files("../bed_files/Tsankov_scATAC/migrated_to_hg19/", 
                                        pattern="IC.*bgz$")
    files = c(files_Tsankov_distal, files_Tsankov_proximal)
    mclapply(files, 
             create_count_overlaps_files,
             metadata=metadata,
             interval_ranges=interval.ranges,
             chain=ch,
             dataset=dataset,
             annotation=annotation,
             which_interval_ranges=which_interval_ranges,
             mc.cores=cores)
  }
  # if (annotation == "Tsankov_separate_fibroblasts") {
  #   if (!file.exists("../metadata/Tsankov_fibro-fibro+C12+fibro+C14.csv")) {
  #     refined_annotation = as_tibble(read.csv("../metadata/metadata_Tsankov_refined_fibroblasts.csv"))
  #     refined_annotation = refined_annotation %>% 
  #                          filter(Clusters %in% c("C12", "C14"))
  #     refined_annotation["Clusters"] = paste("Fibroblasts",
  #                                            refined_annotation[["Clusters"]],
  #                                            sep="_")
  #     idx = match(refined_annotation[["X"]], metadata_tsankov_distal[["X"]])
  #     idx_idx = seq(1, length(idx))
  #     idx_idx = idx_idx[!is.na(idx)]
  #     idx = idx[!is.na(idx)]
  #     metadata_tsankov_distal[idx, "celltypes"] = refined_annotation[idx_idx, 
  #                                                                    "Clusters"]
  #     write.csv(metadata_tsankov_distal, 
  #               "../metadata/Tsankov_fibro-fibro+C12+fibro+C14.csv")
  #   } else {
  #     metadata_tsankov_distal = 
  #       read.csv("../metadata/Tsankov_fibro-fibro+C12+fibro+C14.csv")
  #   }
  # }
  # else if (annotation == "Tsankov_de_novo_Basal_nfrags_filter_1000_tss_filter_4") {
  #   metadata_tsankov_proximal = 
  #     read.csv("../metadata/ArchR_dataset_Tsankov_tissue_all_cell_types_Basal_nfrags_filter_1000_tss_filter_4_annotation.csv")
  # }
  # colnames(metadata_tsankov_proximal)[grepl("celltypes",
  #                                           colnames(metadata_tsankov_proximal))] = "cell_type"
  # colnames(metadata_tsankov_proximal)[grepl("Sample",
  #                                           colnames(metadata_tsankov_proximal))] = "sample"
  # colnames(metadata_tsankov_distal)[grepl("celltypes",
  #                                           colnames(metadata_tsankov_distal))] = "cell_type"
  # colnames(metadata_tsankov_distal)[grepl("Sample",
  #                                           colnames(metadata_tsankov_distal))] = "sample"
  # 
  # 
  # files_Tsankov_distal = list.files("../bed_files/Tsankov_scATAC/migrated_to_hg19/", 
  #                                   pattern="RPL.*bgz$")
  # files_Tsankov_proximal = list.files("../bed_files/Tsankov_scATAC/migrated_to_hg19/", 
  #                                     pattern="IC.*bgz$")
  # dataset_subsets = unlist(strsplit(args$dataset_subsets, ","))
  # if ("proximal" %in% dataset_subsets) {
  #   mclapply(files_Tsankov_proximal, 
  #            create_count_overlaps_files,
  #            metadata=metadata_tsankov_proximal,
  #            interval_ranges=interval.ranges,
  #            chain=ch,
  #            dataset=dataset,
  #            annotation=annotation,
  #            mc.cores=cores)
  # }
  # if ("distal" %in% dataset_subsets) {
  #   mclapply(files_Tsankov_distal,
  #            create_count_overlaps_files,
  #             metadata=metadata_tsankov_distal,
  #             interval_ranges=interval.ranges,
  #             chain=ch,
  #             dataset=dataset,
  #             annotation=annotation,
  #             mc.cores=cores)
  # }
} else if (dataset == "Greenleaf_brain") {
    if (!(file.exists("../metadata/GSE162170_atac_cell_metadata_with_cell_names.txt"))) {
      metadata_greenleaf_brain =
        read.csv("../metadata/GSE162170_atac_cell_metadata.txt.gz",
                 sep="\t")
      cluster_to_cell_names = read.csv("../metadata/scATAC_Colors_greenleaf_brain.txt",
                                       sep = "\t") 
      cell_name_idx = match(metadata_greenleaf_brain[["Iterative.LSI.Clusters"]], 
                         cluster_to_cell_names[["cluster"]])
      cell_types = cluster_to_cell_names[cell_name_idx, "class"]
      metadata_greenleaf_brain["cell_type"] = cell_types
      write.csv(metadata_greenleaf_brain,
                "../metadata/GSE162170_atac_cell_metadata_with_cell_names.txt")
    }
    else {
      if (annotation == "Greenleaf_brain_lowest_level_annotation") {
        metadata_greenleaf_brain =
          read.csv("../metadata/GSE162170_atac_cell_metadata.txt.gz",
                   sep="\t")
        colnames(metadata_greenleaf_brain)[grepl("Iterative.LSI.Clusters", 
                              colnames(metadata_greenleaf_brain))] = "cell_type"
        dir.create("../processed_data/count_overlap_data/Greenleaf_brain_lowest_level_annotation")
      } 
      else {
          metadata_greenleaf_brain =
            read.csv("../metadata/GSE162170_atac_cell_metadata_with_cell_names.txt")
      }
    }
  
  colnames(metadata_greenleaf_brain)[grepl("Sample.ID", colnames(metadata_greenleaf_brain))] = "sample"
  files_greenleaf_brain = list.files("../bed_files/greenleaf_brain_scATAC/migrated_to_hg19", 
                                     pattern=".*fragments\\.tsv\\.bgz$")
  # mclapply(files_greenleaf_brain, create_count_overlaps_files,
  #          cell_number_filter=cell_number_filter,
  #          metadata=metadata_greenleaf_brain,
  #          interval_ranges=interval.ranges,
  #          chain=ch,
  #          dataset=dataset,
  #          mc.cores=cores)
  mclapply(files_greenleaf_brain, create_count_overlaps_files,
           metadata=metadata_greenleaf_brain,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           mc.cores=cores)
} else if (dataset == "Greenleaf_pbmc_bm") {
    if (!(file.exists("../metadata/greenleaf_pbmc_bm.txt"))) {
      metadata_greenleaf_pbmc_bm = readRDS("../metadata/scATAC-Healthy-Hematopoiesis-191120.rds")
      metadata_greenleaf_pbmc_bm = data.frame(barcode = metadata_greenleaf_pbmc_bm$Barcode, 
                                              cell_type = unlist(lapply(
                                              strsplit(metadata_greenleaf_pbmc_bm$BioClassification, 
                                                       split="_"), "[", 2)),
                                              sample = metadata_greenleaf_pbmc_bm$Group)
      
      write.csv(metadata_greenleaf_pbmc_bm, 
                "../metadata/greenleaf_pbmc_bm.txt")
    }
    else {
      metadata_greenleaf_pbmc_bm = 
        read.csv("../metadata/greenleaf_pbmc_bm.txt")
    }
    files_greenleaf_pbmc_bm = list.files("../bed_files/greenleaf_pbmc_bm_scATAC", 
                                         pattern=".*fragments\\.tsv\\.gz")
    mclapply(files_greenleaf_pbmc_bm, 
             create_count_overlaps_files,
             metadata=metadata_greenleaf_pbmc_bm,
             interval_ranges=interval.ranges,
             chain=ch,
             dataset=dataset,
             annotation=annotation,
             which_interval_ranges=which_interval_ranges,
             mc.cores=cores)
} else if (dataset == "Yang_kidney") {
  metadata_Yang = read_excel("../metadata/41467_2021_27660_MOESM4_ESM.xlsx", 
                             skip=1)
  if (which_interval_ranges == "polak") {
    files_Yang = list.files("../bed_files/yang_kidney_scATAC/",
                                pattern = ".*fragments\\.tsv\\.gz")
  }
  else {
    files_Yang = list.files("../bed_files/yang_kidney_scATAC/migrated_to_hg38",
                            pattern = ".*fragments\\.tsv\\.bgz$")
  }
  print(files_Yang)
  colnames(metadata_Yang)[2] = "cell_type"
  colnames(metadata_Yang)[3] = "sample"
  metadata_Yang[, "sample"] = as.character(metadata_Yang$sample)
  metadata_Yang[metadata_Yang[3] == 1, "sample"] = "SRR13679156"
  metadata_Yang[metadata_Yang[3] == 2, "sample"] = "SRR13679157"
  
  mclapply(files_Yang, create_count_overlaps_files,
           metadata=metadata_Yang,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           which_interval_ranges=which_interval_ranges,
           mc.cores=cores)
} else if (dataset == "Greenleaf_colon") {
  if (annotation == "default_annotation") {
    # TODO: MUST FIX NAMING OF SAMPLES
    metadata = read.csv("../metadata/greenleaf_colon_metadata.csv")
  }
  files_colon_greenleaf = list.files("../bed_files/greenleaf_colon_scATAC/migrated_to_hg19",
                                     pattern = ".*fragments\\.tsv\\.bgz$")
  colnames(metadata)[4] = "cell_type"
  colnames(metadata)[2] = "sample"
  # metadata_Yang[, "sample"] = as.character(metadata_Yang$sample)
  # metadata_Yang[metadata_Yang[3] == 1, "sample"] = "SRR13679156"
  # metadata_Yang[metadata_Yang[3] == 2, "sample"] = "SRR13679157"
  
  mclapply(files_colon_greenleaf, create_count_overlaps_files,
           metadata=metadata,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           which_interval_ranges=which_interval_ranges,
           mc.cores=cores)
} else if (dataset == "Rawlins_fetal_lung") {
  if (annotation == "default_annotation") {
    metadata = read.csv("../metadata/rawlins_fetal_lung_metadata.csv")
  }
  files_rawlins_fetal_lung = list.files("../bed_files/rawlins_fetal_lung_scATAC/migrated_to_hg19",
                                     pattern = "\\.tsv\\.bgz$")
  mclapply(files_rawlins_fetal_lung, create_count_overlaps_files,
           metadata=metadata,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           which_interval_ranges=which_interval_ranges,
           mc.cores=cores)
} else if (dataset == "Wang_lung") {
  if (annotation == "default_annotation") {
    if (file.exists("../metadata/wang_adult_metadata.csv")) {
      metadata = read.csv("../metadata/wang_adult_metadata.csv", row.names = 1)
    }
    else {
      metadata = read.table("../metadata/GSE161381_lung_snATAC.UMAP.cluster_labels.txt.gz",
                            sep="\t", header=1)
      temp = strsplit(metadata[["cell"]], split="_")
      metadata["sample"] = unlist(lapply(temp, "[", 1))
      # metadata["cell"] = unlist(lapply(temp, "[", 2))
      metadata = metadata[, c(1, 4, 5)]
      colnames(metadata)[2] = "cell_type"
      write.csv(metadata, "../metadata/wang_adult_metadata.csv")
    }
  }
  files_wang_lung = list.files("../bed_files/wang_adult_lung",
                                        pattern = "D.*")
  mclapply(files_wang_lung, create_count_overlaps_files,
           metadata=metadata,
           interval_ranges=interval.ranges,
           chain=ch,
           dataset=dataset,
           annotation=annotation,
           which_interval_ranges=which_interval_ranges,
           mc.cores=cores)
}
