library(stringr)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)

library(exomeCopy)
mut = readRDS('raw_dir/mutation_data/mutations.aggregated.PCAWG.RK.20181215.Rds')
# get all mutations, without distinctions e.g clonal vs subclonal 
mut = mut[, 1:37]
load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
interval.ranges

# args = commandArgs(trailingOnly=TRUE)
# if (length(args) != 1) {
#   stop("Cell filter argument must be supplied", call.=FALSE)
# } else {
#   CELL_NUMBER_FILTER = args[1]
# } 

dir.create("figures") 

get_per_cancer_mut_data <- function(all_mutations, interval_ranges) {
  cancer_types_mut_data <- list()
  for (cancer_type in colnames(all_mutations)) {
    idx_cancer_type = match(rownames(as.data.frame(interval_ranges)), 
                            rownames(all_mutations[cancer_type]))
    cancer_types_mut_data = append(cancer_types_mut_data,
                                   list(mut[cancer_type][idx_cancer_type, ]))
  }
  return(cancer_types_mut_data)
}

get_tissue_name <- function(filename) {
  if (grepl("frontal_cortex", filename)) {
    tissue_name = "frontal_cortex"
  }
  else {
    tissue_name = unlist(strsplit(str_extract(filename, 
                                              "count_overlaps_.*[0-9]_.*_SM"), 
                                  "_"))
    tissue_name = str_to_title(paste(tissue_name[4:(length(tissue_name)-1)], 
                                     collapse=" "))
  }
  return(tissue_name)
}

save_heatmap_file <- function(path, df, fragment_per_row) {
  pdf(path)
  row_ha = rowAnnotation(fragment_number = anno_text(fragment_per_row))
  h = Heatmap(df, row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 6), 
              left_annotation = row_ha,
              show_column_dend = FALSE, show_row_dend = FALSE)
  print(h)
  dev.off()
}

create_per_tissue_heatmap <- function(cor_df, n, i, frag_df) {
  tissue = n[i]
  tissue = gsub(" ", "_", tolower(tissue), )
  heatmap_filename = tissue
  correlations = cor_df[[i]]
  fragments_per_row = frag_df[[i]]
  
  path = paste("figures/per_tissue_heatmaps/", heatmap_filename, ".pdf", sep="")
  save_heatmap_file(path, correlations, fragments_per_row) 
  
  correlations_rowwise = t(scale(t(correlations)))
  correlations_columnwise = scale(correlations)
  row_heatmap_filename = paste(heatmap_filename, "row_normalized", sep="_")
  column_heatmap_filename = paste(heatmap_filename, "column_normalized", 
                                  sep="_")
  
  row_path = paste("figures/per_tissue_heatmaps/", row_heatmap_filename, ".pdf",
                   sep="")
  column_path = paste("figures/per_tissue_heatmaps/", column_heatmap_filename, 
                      ".pdf", sep="")
  
  save_heatmap_file(row_path, correlations_rowwise, fragments_per_row)
  save_heatmap_file(column_path, correlations_columnwise, fragments_per_row)
}

save_combined_overlaps <- function(filepaths,
                                   combined_filepath,
                                   unsquashed_overlaps_filepath=NULL,
                                   unsquashed_tissue_names_filepath=NULL) {
  combined_count_overlaps = data.frame()
  unsquashed_count_overlaps = data.frame()
  unsquashed_tissue_names = c()
  for (f in filepaths) {
    print(paste("Processing count overlaps for file: ", f, sep=""))
    if (grepl("frontal_cortex", f)) {
      next
    }
    tissue_name <- get_tissue_name(f)
    count_overlaps = readRDS(f)
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
                                   row.names = paste(tissue_name,
                                                     names(count_overlaps)))
    if (grepl("subdivided", f)) {
      colnames(count_overlaps) = rep(rownames(mut), each=5)
    }
    
    if (!is.null(unsquashed_overlaps_filepath)) {
      unsquashed_tissue_names = append(unsquashed_tissue_names, rep(tissue_name,
                                                        dim(count_overlaps)[1]))
      unsquashed_count_overlaps = rbind(unsquashed_count_overlaps, 
                                        count_overlaps)
    }
    
    for (i in 1:nrow(count_overlaps)) {
      count_overlaps_exists = rownames(count_overlaps)[i] %in%
        rownames(combined_count_overlaps)
      tryCatch(
        if (count_overlaps_exists) {
          idx = match(rownames(count_overlaps)[i],
                      rownames(combined_count_overlaps))
          combined_count_overlaps[idx, ] = combined_count_overlaps[idx, ] +
            count_overlaps[i, ]
        }
        else {
          combined_count_overlaps = rbind(combined_count_overlaps,
                                          count_overlaps[i, ])
        },
        error = function(err) {
          print(paste(f, "has no count overlaps"), sep=" ")
        })
    }
  }
  saveRDS(combined_count_overlaps, combined_filepath)
  if (!is.null(unsquashed_overlaps_filepath)) {
    saveRDS(unsquashed_count_overlaps, unsquashed_overlaps_filepath)
    saveRDS(unsquashed_tissue_names, unsquashed_tissue_names_filepath)
  }
}

create_and_save_mutation_df_all_cancers <- function() {
  cancer_types_mut_data = get_per_cancer_mut_data(mut, interval.ranges)
  mut_count_data = as.data.frame(do.call(cbind, cancer_types_mut_data))
  colnames(mut_count_data) = colnames(mut)
  rownames(mut_count_data) = names(interval.ranges)
  write.csv(mut_count_data, "processed_data/mut_count_data.csv")
  return(mut_count_data)
}

mut_count_data = create_and_save_mutation_df_all_cancers()

combined_data_path = "processed_data/count_overlap_data/combined_count_overlaps/"
dir.create(combined_data_path)
combined_filepath = paste(combined_data_path, 
                          "count_filter_",
                          CELL_NUMBER_FILTER, 
                          "_combined_count_overlaps.rds", sep="")

unsquashed_overlaps_filepath = paste(combined_data_path, 
                                   "count_filter_",
                                    CELL_NUMBER_FILTER, 
                                    "_unsquashed_count_overlaps.rds", sep="")

unsquashed_tissue_names_filepath = paste(combined_data_path, 
                                     "count_filter_",
                                     CELL_NUMBER_FILTER, 
                                     "_unsquashed_tissue_names.rds", sep="")

combined_subdivided_filepath = paste(combined_data_path, 
                          "count_filter_",
                          CELL_NUMBER_FILTER, 
                          "_subdivided_combined_count_overlaps.rds", sep="")

#### One heatmap, all tissues ####
if (!file.exists(combined_filepath) || 
    !file.exists(unsquashed_overlaps_filepath)) {
  pattern = paste("count_filter_", CELL_NUMBER_FILTER, "_count_overlaps", 
                  sep="")
  filepaths = list.files("processed_data/count_overlap_data", pattern = pattern, 
                         full.names = TRUE)
  save_combined_overlaps(filepaths, combined_filepath, 
                         unsquashed_overlaps_filepath,
                         unsquashed_tissue_names_filepath)
}

if (!file.exists(combined_subdivided_filepath)) {
  pattern_subdivided = paste("count_filter_", CELL_NUMBER_FILTER, "_subdivided", 
                             sep="")
  filepaths_subdivided = list.files("processed_data/count_overlap_data", 
                                    pattern = pattern_subdivided, 
                                    full.names = TRUE)
  save_combined_overlaps(filepaths_subdivided, combined_subdivided_filepath)
}

#combined_count_overlaps = t(readRDS(combined_subdivided_filepath))

combined_count_overlaps = t(readRDS(combined_filepath))
correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
correlations = t(scale(t(correlations)))

#### Tissue specific heatmaps ####
tissue_specific_counts_path = "processed_data/count_overlap_data/combined_count_overlaps/tissue_specific_counts_list.rds"
tissue_specific_counts = list()

if (!file.exists(tissue_specific_counts_path)) {
  for (f in list.files("processed_data/count_overlap_data/", pattern=pattern,
                          full.names = TRUE)) {
    count_overlaps = readRDS(f)
    tissue_name <- get_tissue_name(f)
    
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
                                   row.names = paste(tissue_name,
                                                     names(count_overlaps)))
    if (dim(count_overlaps)[1] != 0) {
      if (!tissue_name %in% names(tissue_specific_counts)) {
        tissue_specific_counts[[tissue_name]] = count_overlaps
      }
      else {
        for (row in 1:nrow(count_overlaps)) {
          current_cell_subtype = rownames(count_overlaps)[row]
          current_tissue_counts = tissue_specific_counts[[tissue_name]]
          if(current_cell_subtype %in% 
             rownames(current_tissue_counts)) {
            idx = match(current_cell_subtype, rownames(current_tissue_counts)) 
            tissue_specific_counts[[tissue_name]][idx, ] = 
              current_tissue_counts[idx, ] + count_overlaps[row, ]
          }
          else {
            tissue_specific_counts[[tissue_name]][current_cell_subtype, ] = 
              count_overlaps[row, ]
          }
        }
      }
    }
  }
  saveRDS(tissue_specific_counts, tissue_specific_counts_path)
}

tissue_specific_counts = readRDS("processed_data/count_overlap_data/combined_count_overlaps/tissue_specific_counts_list.rds")
tissue_specific_counts = lapply(tissue_specific_counts, t)
per_tissue_cor = lapply(tissue_specific_counts, cor, mut_count_data, 
                        use="complete.obs")
per_tissue_fragment_counts = lapply(tissue_specific_counts, colSums)

dir.create("figures/per_tissue_heatmaps")
lapply(seq_along(per_tissue_cor), create_per_tissue_heatmap, 
       cor_df=per_tissue_cor, 
       n=names(per_tissue_cor),
       frag_df=per_tissue_fragment_counts)


#### Avg correlation vs number of mutations ####
correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
muts_per_cancer_type = apply(mut, 2, sum)
avg_corrs_across_cell_types = apply(correlations, 2, mean)
df = data.frame(muts_per_cancer_type, avg_corrs_across_cell_types)
ggplot(df, aes(x=muts_per_cancer_type, y=avg_corrs_across_cell_types)) + 
  xlab("Number of mutations") +
  ylab("Average Correlation") +
  geom_point() + 
  stat_cor(method = "pearson", label.x = 5000000, label.y = 1.5)

cor(muts_per_cancer_type, avg_corrs_across_cell_types)

#### Min correlation vs number of fragments ####
cell_type = 'Fibroblast'
total_fragments_per_sample = rowSums(unsquashed_count_overlaps)
unsquashed_correlations = cor(t(unsquashed_count_overlaps), mut_count_data, 
                              use = "complete.obs")
min_corr_per_sample = apply(unsquashed_correlations, 1, min)
df = data.frame(fragments=total_fragments_per_sample, 
                correlation=min_corr_per_sample,
                tissue=unsquashed_tissue_names)
df = df[grepl(cell_type, rownames(df), ignore.case = T), ]
ggplot(df, aes(x=fragments, y=correlation, color=tissue)) +
  geom_point() 
# stat_cor(method = "pearson",label.x = 5000000, label.y = -1.5)

cor(df[,'fragments'], df[, 'correlation'])

#### Polak Count overlaps ####
polak_combined_count_overlaps = data.frame()
if (!file.exists("polak_count_overlap_data/processed_count_overlaps/polak_combined_count_overlaps.rds")) {
  polak_combined_count_overlaps = data.frame()
  for (f in list.files("polak_count_overlap_data", full.names=TRUE, pattern=
                          "count_overlaps_")) {
    if (grepl("uniques", f)) {
      next
    }
    else {
      tissue_or_cell_name = unlist(strsplit(f, ".", fixed = TRUE))[2]
      
      if (grepl("Left", tissue_or_cell_name) || 
          grepl("Right", tissue_or_cell_name)) {
        tissue_or_cell_name = as.list(unlist(strsplit(tissue_or_cell_name, "_")))
        tissue_or_cell_name = tissue_or_cell_name[1:length(tissue_or_cell_name)-1]
        tissue_or_cell_name = str_to_title(paste(unlist(tissue_or_cell_name),
                           collapse=" "))
      }
      else {
        tissue_or_cell_name = str_to_title(paste(unlist(strsplit(tissue_or_cell_name, "_")),
                                       collapse=" "))
      }
    }
    count_overlaps = readRDS(f)
    count_overlaps = t(as.data.frame(count_overlaps))
    row.names(count_overlaps) = tissue_or_cell_name
    count_overlaps_exists = rownames(count_overlaps) %in%
                            rownames(polak_combined_count_overlaps)
    tryCatch(
          if (count_overlaps_exists) {
            idx = match(rownames(count_overlaps), rownames(polak_combined_count_overlaps))
            polak_combined_count_overlaps[idx, ] = 
              polak_combined_count_overlaps[idx, ] + count_overlaps[1, ]
          }
          else {
            polak_combined_count_overlaps = rbind(polak_combined_count_overlaps, 
                                                  count_overlaps)
          },
          error = function(err) {
            print(paste(f, "has no count overlaps"), sep=" ")
     })
  }
  saveRDS(polak_combined_count_overlaps, 
          "polak_count_overlap_data/processed_count_overlaps/polak_combined_count_overlaps.rds")
}

polak_combined_filepath = 
  "polak_count_overlap_data/processed_count_overlaps/polak_combined_count_overlaps.rds"

polak_combined_count_overlaps = t(readRDS(polak_combined_filepath))

correlations = cor(polak_combined_count_overlaps, mut_count_data, 
                   use="complete.obs")
correlations = scale(correlations)

pdf ("/home/mdanb/research/mount_sinai/bing_ren/corrs.pdf", width = 10, 
     height = 30)
h = Heatmap(correlations,
           row_names_gp = gpar(fontsize = 6),
           column_names_gp = gpar(fontsize = 6),height=18,
)
print(h)
dev.off()
