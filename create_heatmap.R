library(stringr)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)

mut = readRDS('raw_dir/mutation_data/mutations.aggregated.PCAWG.RK.20181215.Rds')
load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
interval.ranges

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Cell filter argument must be supplied", call.=FALSE)
} else {
  CELL_NUMBER_FILTER = args[1]
} 

dir.create("figures") 

get_tissue_name <- function(filename) {
  if (grepl("frontal_cortex", file)) {
    tissue_name = "frontal_cortex"
  }
  else {
    tissue_name = unlist(strsplit(str_extract(file, 
                                              "count_overlaps_.*[0-9]_.*_SM"), 
                                  "_"))
    tissue_name = str_to_title(paste(tissue_name[4:(length(tissue_name)-1)], 
                                     collapse=" "))
  }
  return(tissue_name)
}

save_heatmap_file <- function(path, df) {
  pdf(path)
  h = Heatmap(df, row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 6))
  print(h)
  dev.off()
}

create_per_tissue_heatmap <- function(cor_df, n, i) {
  tissue = n[i]
  tissue = gsub(" ", "_", tolower(tissue), )
  heatmap_filename = tissue
  correlations = cor_df[[i]]
  
  path = paste("figures/per_tissue_heatmaps/", heatmap_filename, ".pdf", sep="")
  save_heatmap_file(path, correlations) 
  
  correlations_rowwise = t(scale(t(correlations)))
  correlations_columnwise = scale(correlations)
  row_heatmap_filename = paste(heatmap_filename, "row_normalized", sep="_")
  column_heatmap_filename = paste(heatmap_filename, "column_normalized", 
                                  sep="_")
  
  row_path = paste("figures/per_tissue_heatmaps/", row_heatmap_filename, ".pdf",
                   sep="")
  column_path = paste("figures/per_tissue_heatmaps/", column_heatmap_filename, 
                      ".pdf", sep="")
  
  save_heatmap_file(row_path, correlations_rowwise)
  save_heatmap_file(column_path, correlations_columnwise)
}
# get all mutations, without distinctions e.g clonal vs subclonal 
mut = mut[, 1:37]

cancer_types_mut_data = list()

for (cancer_type in colnames(mut)) {
  idx_cancer_type = match(rownames(as.data.frame(interval.ranges)), 
                          rownames(mut[cancer_type]))
  cancer_types_mut_data = append(cancer_types_mut_data,
                                 list(mut[cancer_type][idx_cancer_type, ]))
}

mut_count_data = as.data.frame(do.call(cbind, cancer_types_mut_data))
colnames(mut_count_data) = colnames(mut)
rownames(mut_count_data) = names(interval.ranges)

combined_filepath = paste("count_overlap_data/processed_count_overlaps/", 
                          "count_filter_",
                          CELL_NUMBER_FILTER, 
                          "_combined_count_overlaps.rds", sep="")

unsquashed_overlaps_filepath = paste("count_overlap_data/processed_count_overlaps/", 
                                   "count_filter_",
                                    CELL_NUMBER_FILTER, 
                                    "_unsquashed_count_overlaps.rds", sep="")

unsquashed_tissue_names_filepath = paste("count_overlap_data/processed_count_overlaps/", 
                                     "count_filter_",
                                     CELL_NUMBER_FILTER, 
                                     "_unsquashed_tissue_names.rds", sep="")

pattern = paste("count_filter_", CELL_NUMBER_FILTER, "_", sep="")
if (!file.exists(combined_filepath) || !file.exists(unsquashed_overlaps_filepath)) {
  combined_count_overlaps = data.frame()
  unsquashed_count_overlaps = data.frame()
  unsquashed_tissue_names = c()
  for (file in list.files("count_overlap_data", pattern = pattern, 
                          full.names = TRUE)) {
    tissue_name <- get_tissue_name(file)
    count_overlaps = readRDS(file)
    count_overlaps = as.data.frame(do.call(rbind, count_overlaps), 
                                   row.names = paste(tissue_name, 
                                                     names(count_overlaps)))
    unsquashed_tissue_names = append(unsquashed_tissue_names, rep(tissue_name, 
                                                                  dim(count_overlaps)[1]))
    unsquashed_count_overlaps = rbind(unsquashed_count_overlaps, count_overlaps)
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
          print(paste(file, "has no count overlaps"), sep=" ")
        })
    }
  }
  saveRDS(combined_count_overlaps, combined_filepath)
  saveRDS(unsquashed_count_overlaps, unsquashed_overlaps_filepath)
  saveRDS(unsquashed_tissue_names, unsquashed_tissue_names_filepath)
}

combined_filepath = paste("count_overlap_data/processed_count_overlaps/", 
                          "count_filter_",
                          100, 
                          "_combined_count_overlaps.rds", sep="")

combined_count_overlaps = t(readRDS(combined_filepath))

correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
correlations = t(scale(t(correlations)))
# correlations = scale(correlations)

#write.csv(correlations, "all_corrs.csv")

#png(file="/home/mdanb/research/mount_sinai/bing_ren/corrs.png")
#pdf ("/home/mdanb/research/mount_sinai/bing_ren/corrs.pdf", width = 10, height = 30)
#h = Heatmap(correlations,
#            row_names_gp = gpar(fontsize = 6),
#            column_names_gp = gpar(fontsize = 6),height=18,
#)
#print(h)
#dev.off()


tissue_specific_counts_path = "count_overlap_data/processed_count_overlaps/tissue_specific_counts_list.rds"
tissue_specific_counts = list()

if (!file.exists(tissue_specific_counts_path)) {
  for (file in list.files("count_overlap_data", pattern=pattern,
                          full.names = TRUE)) {
    count_overlaps = readRDS(file)
    tissue_name <- get_tissue_name(file)
    
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
            tissue_specific_counts[[tissue_name]][idx, ] = current_tissue_counts[idx, ] +
              count_overlaps[row, ]
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

tissue_specific_counts = readRDS("count_overlap_data/processed_count_overlaps/tissue_specific_counts_list.rds")
tissue_specific_counts = lapply(tissue_specific_counts, t)
per_tissue_cor = lapply(tissue_specific_counts, cor, mut_count_data, 
                        use="complete.obs")

dir.create("figures/per_tissue_heatmaps")
lapply(seq_along(per_tissue_cor), create_per_tissue_heatmap, 
       cor_df=per_tissue_cor, 
       n=names(per_tissue_cor))


#### ####
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

#### ####
# obtain_x_y_group_triplets <- function(df, n, i, cell_type, 
#                                       combined_count_overlaps) {
#   df = df[[i]]
#   tissue = n[i]
#   pattern = paste(tissue, ".*", cell_type, sep="")
#   x = sum(combined_count_overlaps[, grepl(pattern=pattern, 
#                                   colnames(combined_count_overlaps))])
#   relevant_cells = df[grepl(cell_type, rownames(df), ignore.case = T), ]
#   if (!is.null(dim(relevant_cells))) {
#     y = apply(relevant_cells, 1,  min)
#   }
#   else {
#     y = min(relevant_cells)
#   }
#   return(data.frame(cell_num=rep(x, length(y)), corr=y, tissue=tissue,
#                     row.names=NULL))
# }

# fibroblast_min_corrs_per_tissue = 
#                    lapply(seq_along(per_tissue_cor), obtain_x_y_group_triplets,
#                    df = per_tissue_cor, n = names(per_tissue_cor), 
#                    cell_type ='Fibroblast', 
#                    combined_count_overlaps = combined_count_overlaps) %>%
#                    bind_rows()

# cor(fibroblast_min_corrs_per_tissue['cell_num'], fibroblast_min_corrs_per_tissue['corr'])
cell_type = 'Fibroblast'
total_fragments_per_sample = rowSums(unsquashed_count_overlaps)
unsquashed_correlations = cor(t(unsquashed_count_overlaps), mut_count_data, 
                              use = "complete.obs")
min_corr_per_sample = apply(unsquashed_correlations, 1, min)
df = data.frame(fragments=total_fragments_per_sample, correlation=min_corr_per_sample,
                tissue=unsquashed_tissue_names)
df = df[grepl(cell_type, rownames(df), ignore.case = T), ]
ggplot(df, aes(x=fragments, y=correlation, color=tissue)) +
  geom_point() 

cor(df[,'fragments'], df[, 'correlation'])

#+
 # stat_cor(method = "pearson",label.x = 5000000, label.y = -1.5)

#### ####


# polak_combined_count_overlaps = data.frame()
# if (!file.exists("polak_count_overlap_data/polak_combined_count_overlaps.rds")) {
#   polak_combined_count_overlaps = data.frame()
#   for (file in list.files("polak_count_overlap_data", full.names=TRUE)) {
#     if (grepl("IMR90", file)) {
#       tissue_name = "IMR90"
#     }
#     else {
#       tissue_name = unlist(strsplit(file, ".", fixed = TRUE))[2]
#       tissue_name = str_to_title(paste(unlist(strsplit(mystr, "_")), 
#                                        collapse=" "))
#     }
#     count_overlaps = readRDS(file)
#     count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
#                                    row.names = paste(tissue_name,
#                                                      names(count_overlaps)))
#     for (i in 1:nrow(count_overlaps)) {
#       count_overlaps_exists = rownames(count_overlaps)[i] %in%
#         rownames(combined_count_overlaps)
#       tryCatch(
#         if (count_overlaps_exists) {
#           idx = match(rownames(count_overlaps)[i], rownames(combined_count_overlaps))
#           combined_count_overlaps[idx, ] = combined_count_overlaps[idx, ] +
#             count_overlaps[i, ]
#         }
#         else {
#           combined_count_overlaps = rbind(combined_count_overlaps, count_overlaps[i, ])
#         },
#         error = function(err) {
#           print(paste(file, "has no count overlaps"), sep=" ")
#         })
#     }
#   }
#   saveRDS(combined_count_overlaps, "count_overlap_data/combined_count_overlaps.rds")
# }
