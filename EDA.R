library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)
library(stringr)
library(exomeCopy)
library(parallel)
library(RColorBrewer)
library(gtools)
source("utils.R")

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
dir.create("figures") 
CELL_NUMBER_FILTER = 1

load_mutation_data <- function() {
  mut = readRDS('raw_dir/mutation_data/mutations.aggregated.PCAWG.RK.20181215.Rds')
  # get all mutations, without distinctions e.g clonal vs subclonal 
  mut = mut[, 1:37]
  return(mut)
}

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

save_heatmap_file <- function(path, width, height, df, ...) {
  pdf(path, width=width, height=height)
  h = Heatmap(df, show_column_dend = FALSE, show_row_dend = FALSE, ...)
  print(h)
  dev.off()
}

create_and_save_mutation_df_all_cancers <- function() {
  if (!file.exists("processed_data/mut_count_data.csv")) {
    cancer_types_mut_data = get_per_cancer_mut_data(mut, interval.ranges)
    mut_count_data = as.data.frame(do.call(cbind, cancer_types_mut_data))
    colnames(mut_count_data) = colnames(mut)
    rownames(mut_count_data) = names(interval.ranges)
    write.csv(mut_count_data, "processed_data/mut_count_data.csv")
  }
  else {
    mut_count_data = read.csv("processed_data/mut_count_data.csv")
  }
  return(mut_count_data)
}

mut = load_mutation_data()
mut_count_data = create_and_save_mutation_df_all_cancers()
rownames(mut_count_data) = mut_count_data[, 1]
mut_count_data = mut_count_data[, 2:length(colnames(mut_count_data))]
  
combined_data_path = "processed_data/count_overlap_data/combined_count_overlaps/"

combined_filepath = paste(combined_data_path,
                          "count_filter_",
                          CELL_NUMBER_FILTER,
                          "_combined_count_overlaps.rds", sep="")
combined_count_overlaps = t(readRDS(combined_filepath))

combined_filepath_shendure = paste(combined_data_path,
                                   "shendure_count_filter_",
                                   CELL_NUMBER_FILTER,
                                   "_combined_count_overlaps.rds", sep="")

combined_filepath_tsankov = paste(combined_data_path,
                                   "tsankov_count_filter_",
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

# #### Cell type Cell type correlation ####
# # Expect cells of same type from different tissues to cluster together, rather
# # than cells from same tissue to cluster 
# cell_type_cell_type_corr = cor(combined_count_overlaps, combined_count_overlaps)
# total_fragments = apply(combined_count_overlaps, 2, sum)
# ha = rowAnnotation(fragment_counts = anno_barplot(log2(total_fragments)))
# save_heatmap_file("figures/cell_type_cell_type_corr.pdf",
#                   30, 30, cell_type_cell_type_corr, 
#                   row_names_gp = gpar(fontsize = 6),
#                   column_names_gp = gpar(fontsize = 6), left_annotation = ha)
# 
# save_heatmap_file("figures/cell_type_cell_type_corr_pearson.pdf",
#                   30, 30, cell_type_cell_type_corr, 
#                   row_names_gp = gpar(fontsize = 6),
#                   column_names_gp = gpar(fontsize = 6),
#                   left_annotation = ha, 
#                   clustering_distance_rows = "pearson",
#                   clustering_distance_columns = "pearson")
# 
# #### One heatmap, all tissues ####
# correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
# correlations = t(scale(t(correlations)))
# correlations = scale(correlations)
# save_heatmap_file("figures/all_cells_corr_double_normalized.pdf",
#                   30, 30, correlations, row_names_gp = gpar(fontsize = 6),
#                   column_names_gp = gpar(fontsize = 6))
# 
# #### Tissue specific heatmaps ####
# combine_count_overlaps_per_tissue <- function(f, tissue_specific_counts) {
#   count_overlaps = readRDS(f)
#   tissue_name <- get_tissue_name(f)
#   
#   count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
#                                  row.names = paste(tissue_name,
#                                                    names(count_overlaps)))
#   if (dim(count_overlaps)[1] != 0) {
#     if (!tissue_name %in% names(tissue_specific_counts)) {
#       tissue_specific_counts[[tissue_name]] = count_overlaps
#     }
#     else {
#       for (row in 1:nrow(count_overlaps)) {
#         current_cell_subtype = rownames(count_overlaps)[row]
#         current_tissue_counts = tissue_specific_counts[[tissue_name]]
#         if(current_cell_subtype %in% 
#            rownames(current_tissue_counts)) {
#           idx = match(current_cell_subtype, rownames(current_tissue_counts)) 
#           tissue_specific_counts[[tissue_name]][idx, ] = 
#             current_tissue_counts[idx, ] + count_overlaps[row, ]
#         }
#         else {
#           tissue_specific_counts[[tissue_name]][current_cell_subtype, ] = 
#             count_overlaps[row, ]
#         }
#       }
#     }
#   }
#   return(tissue_specific_counts)
# }
# 
# create_per_tissue_heatmap <- function(cor_df, n, i, frag_df) {
#   tissue = n[i]
#   tissue = gsub(" ", "_", tolower(tissue), )
#   heatmap_filename = tissue
#   correlations = cor_df[[i]]
#   fragments_per_row = frag_df[[i]]
#   
#   path = paste("figures/per_tissue_heatmaps/", heatmap_filename, ".pdf", sep="")
#   ha = rowAnnotation(fragment_counts = anno_barplot(fragments_per_row))
#   save_heatmap_file(path, 40, 40, correlations, left_annotation=ha) 
#   
#   correlations_rowwise = t(scale(t(correlations)))
#   correlations_columnwise = scale(correlations)
#   row_heatmap_filename = paste(heatmap_filename, "row_normalized", sep="_")
#   column_heatmap_filename = paste(heatmap_filename, "column_normalized", 
#                                   sep="_")
#   
#   row_path = paste("figures/per_tissue_heatmaps/", row_heatmap_filename, ".pdf",
#                    sep="")
#   column_path = paste("figures/per_tissue_heatmaps/", column_heatmap_filename, 
#                       ".pdf", sep="")
#   
#   save_heatmap_file(row_path, 40, 40, correlations_rowwise, 
#                     left_annotation=ha)
#   save_heatmap_file(column_path, 40, 40, correlations_columnwise, 
#                     left_annotation=ha)
# }
# 
# save_tissue_specific_heatmaps <- function(pattern_of_files_to_consider, 
#                                           tissue_specific_counts_path) {
#   tissue_specific_counts = list()
#   if (!file.exists(tissue_specific_counts_path)) {
#     # pattern = paste("^count_filter", cell_number_filter, sep="_")
#     for (f in list.files("processed_data/count_overlap_data", 
#                          pattern=pattern_of_files_to_consider,
#                          full.names = TRUE)) {
#       tissue_specific_counts <- combine_count_overlaps_per_tissue(f, 
#                                                                   tissue_specific_counts)
#     }
#     saveRDS(tissue_specific_counts, tissue_specific_counts_path)
#   }
# }
# 
# tissue_specific_counts_filename = paste("tissue_specific_counts_list", 
#                                         "count_filter",
#                                         CELL_NUMBER_FILTER,
#                                         sep = "_")
# tissue_specific_counts_filename = paste0(tissue_specific_counts_filename, 
#                                          ".rds") 
# tissue_specific_counts_path = paste("processed_data",
#                                     "count_overlap_data",
#                                     "combined_count_overlaps",
#                                     tissue_specific_counts_filename,
#                                     sep="/")
# 
# pattern_of_files_to_consider = paste("^count_filter", 
#                                      CELL_NUMBER_FILTER, 
#                                      "count_overlaps.*",
#                                      sep="_")
# 
# save_tissue_specific_heatmaps(pattern_of_files_to_consider,
#                               tissue_specific_counts_path)
# 
# # if (!file.exists(tissue_specific_counts_path)) {
# #   pattern = paste("^count_filter", CELL_NUMBER_FILTER, sep="_")
# #   for (f in list.files("processed_data/count_overlap_data", pattern=pattern,
# #                         full.names = TRUE)) {
# #     tissue_specific_counts <- combine_count_overlaps_per_tissue(f, 
# #                                                          tissue_specific_counts)
# #   }
# #   saveRDS(tissue_specific_counts, tissue_specific_counts_path)
# # }
# 
# tissue_specific_counts = readRDS(tissue_specific_counts_path)
# tissue_specific_counts = lapply(tissue_specific_counts, t)
# per_tissue_cor = lapply(tissue_specific_counts, cor, mut_count_data, 
#                         use="complete.obs")
# per_tissue_fragment_counts = lapply(tissue_specific_counts, colSums)
# 
# dir.create("figures/per_tissue_heatmaps")
# lapply(seq_along(per_tissue_cor), create_per_tissue_heatmap, 
#        cor_df=per_tissue_cor, 
#        n=names(per_tissue_cor),
#        frag_df=per_tissue_fragment_counts)
# 
# #### Avg correlation vs number of mutations ####
# correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
# muts_per_cancer_type = apply(mut, 2, sum)
# avg_corrs_across_cell_types = apply(correlations, 2, mean)
# df = data.frame(muts_per_cancer_type, avg_corrs_across_cell_types)
# ggplot(df, aes(x=muts_per_cancer_type, y=avg_corrs_across_cell_types)) + 
#   xlab("Number of mutations") +
#   ylab("Average Correlation") +
#   geom_point() + 
#   stat_cor(method = "pearson", label.x = 5000000, label.y = 1.5)
# 
# #### Min correlation vs number of fragments for Fibroblast ####
# cell_type = 'Fibroblast'
# unsquashed_count_overlaps = readRDS("processed_data/count_overlap_data/combined_count_overlaps/count_filter_100_unsquashed_count_overlaps.rds")
# unsquashed_tissue_names =  readRDS("processed_data/count_overlap_data/combined_count_overlaps/count_filter_100_unsquashed_tissue_names.rds")
# total_fragments_per_sample = rowSums(unsquashed_count_overlaps)
# unsquashed_correlations = cor(t(unsquashed_count_overlaps), mut_count_data, 
#                               use = "complete.obs")
# min_corr_per_sample = apply(unsquashed_correlations, 1, min)
# df = data.frame(fragments=total_fragments_per_sample, 
#                 correlation=min_corr_per_sample,
#                 tissue=unsquashed_tissue_names)
# df = df[grepl(cell_type, rownames(df), ignore.case = T), ]
# corr = round(cor(df[,'fragments'], df[, 'correlation']), 2)
# ggplot(df, aes(x=fragments, y=correlation, color=tissue)) +
#   geom_point() +
#   annotate("text", x=20000000, y=-0.2, label= toString(str_glue("R=", corr)))
#   #stat_cor(method = "pearson",label.x = 5000000, label.y = -1.5)
# 
# #### Polak Count overlaps ####
# polak_combined_count_overlaps = data.frame()
# polak_combined_filepath = 
#   "processed_data/polak_count_overlap_data/combined_count_overlaps/polak_combined_count_overlaps.rds"
# 
# if (!file.exists(polak_combined_filepath)) {
#   polak_combined_count_overlaps = data.frame()
#   for (f in list.files("processed_data/polak_count_overlap_data", 
#                        full.names=TRUE, pattern=
#                           "count_overlaps_")) {
#     if (grepl("uniques", f)) {
#       next
#     }
#     else {
#       tissue_or_cell_name = unlist(strsplit(f, ".", fixed = TRUE))[2]
#       
#       if (grepl("Left", tissue_or_cell_name) || 
#           grepl("Right", tissue_or_cell_name)) {
#         tissue_or_cell_name = as.list(unlist(strsplit(tissue_or_cell_name, "_")))
#         tissue_or_cell_name = tissue_or_cell_name[1:length(tissue_or_cell_name)-1]
#         tissue_or_cell_name = str_to_title(paste(unlist(tissue_or_cell_name),
#                            collapse=" "))
#       }
#       else {
#         tissue_or_cell_name = str_to_title(paste(unlist(strsplit(tissue_or_cell_name, "_")),
#                                        collapse=" "))
#       }
#     }
#     count_overlaps = readRDS(f)
#     count_overlaps = t(as.data.frame(count_overlaps))
#     row.names(count_overlaps) = tissue_or_cell_name
#     count_overlaps_exists = rownames(count_overlaps) %in%
#                             rownames(polak_combined_count_overlaps)
#     tryCatch(
#           if (count_overlaps_exists) {
#             idx = match(rownames(count_overlaps), 
#                         rownames(polak_combined_count_overlaps))
#             polak_combined_count_overlaps[idx, ] = 
#               polak_combined_count_overlaps[idx, ] + count_overlaps[1, ]
#           }
#           else {
#             polak_combined_count_overlaps = rbind(polak_combined_count_overlaps, 
#                                                   count_overlaps)
#           },
#           error = function(err) {
#             print(paste(f, "has no count overlaps"), sep=" ")
#      })
#   }
#   saveRDS(polak_combined_count_overlaps, polak_combined_filepath)
# }
# 
# polak_combined_count_overlaps = t(readRDS(polak_combined_filepath))
# correlations = cor(polak_combined_count_overlaps, mut_count_data, 
#                    use="complete.obs")
# correlations = scale(correlations)
# save_heatmap_file("figures/polak_cell_type_to_cancer_correlations.pdf", 
#                   10, 30, correlations, 
#                   row_names_gp = gpar(fontsize = 6),
#                   column_names_gp = gpar(fontsize = 6))
# 
# #### Our data vs Polak correlation ####
# correlations = t(scale(t(cor(combined_count_overlaps, 
#                              polak_combined_count_overlaps))))
# correlations = scale(correlations)
# save_heatmap_file("figures/our_data_vs_polak_corr.pdf", 10, 30, correlations,
#                   row_names_gp = gpar(fontsize = 6),
#                   column_names_gp = gpar(fontsize = 6))
# 
# #### Cell number / Fragment counts vs R boxplots ####
# correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
# total_fragments = apply(combined_count_overlaps, 2, sum)
# cell_types = rownames(correlations)
# correlations = as_tibble(correlations)
# correlations = correlations %>%
#                add_column(total_fragments=total_fragments) %>%
#                add_column(cell_type=cell_types)
# 
# correlations_long = pivot_longer(correlations, cols = 
#                                    -(total_fragments:cell_type),
#                                  names_to = "cancer_type", 
#                                  values_to = "correlation")
# 
# # correlations_long$total_fragments = 
#   # factor(correlations_long$total_fragments, 
#   #        levels = unique(correlations_long$total_fragments[order (correlations_long$total_fragments)]))
# # geom_boxplot(aes(x=total_fragments, y=correlation)) +
# # theme(axis.text.x=element_text(size=3, angle = 90)) +
# # # xlab("Num Fragments") +
# # scale_x_discrete(name="Num Fragments") +
# 
# medians <- correlations_long %>% 
#   group_by(total_fragments) %>% 
#   summarise(median(correlation))
# 
# colnames(medians)[2] <- "corr_medians"
# correlations_long <- correlations_long %>% right_join(medians)
# 
# ggplot(correlations_long) +
#   geom_boxplot(aes(x=as.factor(total_fragments), y=correlation)) +
#   theme(axis.text.x=element_text(size=3, angle = 90)) +
#   xlab("Num Fragments") 
# 
# correlations_long_no_outliers = correlations_long[correlations_long$corr_medians 
#                                                   < -0.2, ]
# medians = medians[medians$corr_medians < -0.2, ]
# coefs <- coef(lm(corr_medians ~ total_fragments, data = medians))
# 
# ggplot(correlations_long_no_outliers) +
#   geom_boxplot(aes(x=as.factor(total_fragments), y=correlation)) +
#   theme(axis.text.x=element_text(size=3, angle = 90)) +
#   xlab("Num Fragments") +
#   geom_abline(intercept = coefs[1], slope = coefs[2])


#### Fragment counts / TSS vs R for Melanoma boxplots ####
# Fragment counts
run_one_simulation <- function(n_samples, count_overlaps, 
                               sampling_vec, mutations) {
  subsampled_count_overlaps <- rep(0, length(count_overlaps))
  names(subsampled_count_overlaps) <- names(count_overlaps)
  samples <- sample(sampling_vec, n_samples)
  idx = match(names(table(samples)), names(subsampled_count_overlaps))
  subsampled_count_overlaps[idx] = table(samples)
  correlation = cor(subsampled_count_overlaps, mutations, 
                    use="complete.obs")
  return(correlation)
}

create_fragments_subsampling_range <- function(lower, upper, log=T) {
  if (log) {
    lower = log(lower, base=2)
    upper = log(upper, base=2)
  }
  num_fragments_to_subsample <- seq(lower, upper, (upper - lower) / 9)
  # num_fragments_to_subsample[2:(length(num_fragments_to_subsample) - 1)] = 
  #   num_fragments_to_subsample[2:(length(num_fragments_to_subsample) 
  #                                            - 1)]
  if (!log) {
    num_fragments_to_subsample[2:(length(num_fragments_to_subsample) 
                                  - 1)] =
      as.integer(num_fragments_to_subsample[2:(length(num_fragments_to_subsample)
                                               - 1)])
  }
  return(num_fragments_to_subsample)
}

run_simulations_per_num_fragments <- function(num_fragments, log) {
  if (log) {
    num_fragments = as.integer(2^num_fragments)
  }
  samples_10 = mclapply(rep(i, 10), run_one_simulation, 
                       cell_type_count_overlaps,
                       sampling_vec, mutations,
                       mc.cores=8)
  return(unlist(samples_10))
}

run_subsampling_simulations <- function(num_fragments_to_subsample, 
                                        cell_type_count_overlaps, 
                                        sampling_vec, mutations, log) {
  ptm <- proc.time()
  subsampling_simulations <- list()
  list_idx = 1
  set.seed(42)
  
  # subsampling_simulations = mclapply(num_fragments_to_subsample, 
  #                                    run_simulations_per_num_fragments,
  #                                    log, mc.cores=8)
  for (i in num_fragments_to_subsample) {
    if (log) {
      i = as.integer(2^i)
    }
    samples_10 = mclapply(rep(i, 10), run_one_simulation, 
                           cell_type_count_overlaps,
                           sampling_vec, mutations,
                           mc.cores=8)
    subsampling_simulations[[list_idx]] = unlist(samples_10)
    list_idx = list_idx + 1
  }
  names(subsampling_simulations) = num_fragments_to_subsample
  
  if (log) {
    names(subsampling_simulations) = as.integer(unlist(lapply(2, `^`, 
                                                num_fragments_to_subsample)))
  }
  subsampling_simulations = as_tibble(subsampling_simulations)
  print(paste0("time to subsample: ", round((proc.time() - ptm)[["elapsed"]]),
               " seconds"))
  return(subsampling_simulations)
}

add_escape_helper <- function(char_to_escape, cell_type) {
  escaped = paste0("\\", char_to_escape)
  left_idx = unlist(gregexpr(escaped, cell_type))[1]
  left_substr = substr(cell_type, 1, left_idx - 1)
  right_substr = substr(cell_type, left_idx, nchar(cell_type))
  cell_type_for_grep = paste0(left_substr, "\\", right_substr)
  return(cell_type_for_grep)
}

add_escape_if_necessary <- function(cell_type) {
  if (grepl("\\(", cell_type)) {
    cell_type = add_escape_helper("(", cell_type)
  }
  if (grepl("\\)", cell_type)) {
    cell_type = add_escape_helper(")", cell_type)
  }
  if (grepl("\\+", cell_type)) {
    cell_type = add_escape_helper("+", cell_type)
  }
  return(cell_type)
}

plot_and_save_boxplots <- function(correlations_long, cell_types, 
                                   plot_filename, 
                                   correlations_for_tsse_filtered_cells) {
  # colors = get_n_colors(length(cell_types), 4)
  x_for_filtered_correlations = unique(correlations_long %>% 
                                         pull("num_fragments"))
  x_for_filtered_correlations = x_for_filtered_correlations[1:length(correlations_for_tsse_filtered_cells)]
  levels_for_filtered_correlations = x_for_filtered_correlations
  ggplot() +
    geom_point(aes(x=factor(x_for_filtered_correlations,
                            levels=x_for_filtered_correlations),
                   y=correlations_for_tsse_filtered_cells, color="tss filtered"),
               size=3) +
    geom_boxplot(data=correlations_long,
                 aes(x=factor(num_fragments, levels=unique(num_fragments)), 
                     y=correlation, color=cell_type)) +
    theme(axis.text.x=element_text(size=10, angle = 90)) +
    xlab("Num Fragments") +
    ylab("Correlation") +
    labs(color="Cell type") + 
    scale_x_discrete(limits = as.factor(sort(as.integer(
                     unique(correlations_long["num_fragments"]) %>% pull)))) #+
    #scale_color_manual(values=colors)
  # ggplot() +
  #   geom_boxplot(data = subset(correlations_long, cell_type == cell_types[1]), 
  #                aes(x=factor(num_fragments, levels=unique(num_fragments)), 
  #                    y=correlation, color=cell_types[1])) +
  #   geom_boxplot(data = subset(correlations_long, cell_type == cell_types[2]),
  #                aes(x=factor(num_fragments, levels=unique(num_fragments)),
  #                    y=correlation, color=cell_types[2])) +
    
    # geom_boxplot(data = subset(correlations_long, cell_type == cell_types[3]),
    #              aes(x=factor(num_fragments, levels=unique(num_fragments)),
    #                  y=correlation, color=cell_types[3])) +
  ggsave(paste("figures", plot_filename, sep="/"), width = 20, height = 12)
}

remove_bigger_than_curr_fragments <- function(fragments_subsampling_range,
                                              curr_num_fragments,
                                              log) {
  fragments_subsampling_range_reduced <- c()
  for (num_fragments in fragments_subsampling_range) {
    if (log) {
      if (as.integer(2^num_fragments) <= curr_num_fragments) {
        fragments_subsampling_range_reduced <- append(fragments_subsampling_range_reduced,
                                                      num_fragments)
      }
    }
    else {
      if (num_fragments <= curr_num_fragments) {
        fragments_subsampling_range_reduced <- append(fragments_subsampling_range_reduced,
                                                      num_fragments)
      }
    }
  }
  return(fragments_subsampling_range_reduced)
}

get_long_correlations_per_cell_type <- function(i,
                                                cell_types, 
                                                cell_types_total_fragments,
                                                fragments_subsampling_range, 
                                                mutations, 
                                                combined_count_overlaps, 
                                                log) {
  cell_type_total_fragments = cell_types_total_fragments[[i]]
  # cell_type_for_grep = add_escape_if_necessary(cell_types[i])
  # mutations = mut_count_data[cancer_type]
  # cell_type_col_idx = grep(cell_type_for_grep, 
  #                          colnames(combined_count_overlaps), 
  #                          ignore.case=T)
  cell_type_count_overlaps = combined_count_overlaps[, i]
  # cell_type_total_fragments = sum(cell_type_count_overlaps)
  cell_type_fragments_subsampling_range = remove_bigger_than_curr_fragments(fragments_subsampling_range,
                                                                  cell_type_total_fragments,
                                                                  log)
  # upper = cell_type_total_fragments
  sampling_vec = rep(names(cell_type_count_overlaps), 
                     cell_type_count_overlaps)
  # fragments_subsampling_range <- create_fragments_subsampling_range(lower, 
  #                                                                   upper, 
  #                                                                   log)
  subsampling_simulations <- run_subsampling_simulations(cell_type_fragments_subsampling_range, 
                                                         cell_type_count_overlaps, 
                                                         sampling_vec, 
                                                         mutations, 
                                                         log)
  curr_correlations_long = pivot_longer(subsampling_simulations, 
                                        cols = everything(),
                                        names_to = "num_fragments", 
                                        values_to = "correlation")
  curr_correlations_long["cell_type"] = cell_types[i]
  # correlations_long = rbind(correlations_long, curr_correlations_long)
  # num_fragments_per_cell_type = append(num_fragments_per_cell_type,
  #                                      cell_type_total_fragments)
  return(curr_correlations_long)
}

prep_boxplots_per_cancer_type <- function(combined_count_overlaps, 
                                          mut_count_data, cancer_type, 
                                          cell_types, lower, log, 
                                          plot_filename, 
                                          correlations_for_tsse_filtered_cells) {
  correlations_long = tibble()
  cell_type_for_grep = lapply(cell_types, add_escape_if_necessary)
  cell_type_col_idx = lapply(cell_type_for_grep, grep,
                             colnames(combined_count_overlaps), ignore.case=T)
  cell_types_count_overlaps = combined_count_overlaps[, unlist(cell_type_col_idx)]
  cell_types_total_fragments = apply(cell_types_count_overlaps, 2, sum)
  # num_fragments_per_cell_type = apply(combined_count_overlaps, 2, sum)
  # mutations = mut_count_data[cancer_type]
  fragments_subsampling_range <- create_fragments_subsampling_range(lower, 
                                                                    max(cell_types_total_fragments), 
                                                                    log)
  mutations <- mut_count_data[cancer_type]
  cell_types_underscore_sep = gsub(" ", "_", cell_types)
  cell_types_being_considered_combined = paste(cell_types_underscore_sep, 
                                               collapse="_and_")
  correlations_filename = paste(cancer_type, 
                                cell_types_being_considered_combined,
                                sep = "_")
  correlations_filename = paste0(correlations_filename, ".rds")
  correlations_filepath = paste("processed_data", "tsse_filtered_correlations",
                                correlations_filename, sep="/")
  # correlations_filepath = ""
  if (!file.exists(correlations_filepath)) {
    correlations_long = mclapply(seq_along(cell_types),
                                 get_long_correlations_per_cell_type,
                                 colnames(cell_types_count_overlaps),
                                 cell_types_total_fragments,
                                 fragments_subsampling_range,
                                 mutations,
                                 cell_types_count_overlaps,
                                 log,
                                 mc.cores=8)
    correlations_long = bind_rows(correlations_long)
    saveRDS(correlations_long, correlations_filepath)
  }
  else {
    correlations_long = readRDS(correlations_filepath)
  }

  plot_and_save_boxplots(correlations_long, colnames(cell_types_count_overlaps),
                         plot_filename, correlations_for_tsse_filtered_cells)
}

# prep_boxplots <- function(cancer_types) {
#   mclapply(cancer_types, 
#            )
# }

# combine_scATAC <- function(combined_filepath, combined_filepath_shendure,
#                            combined_filepath_tsankov) {
#   combined_count_overlaps = t(readRDS(combined_filepath))
#   shendure_combined_count_overlaps = t(readRDS(combined_filepath_shendure))
#   tsankov_combined_count_overlaps = t(readRDS(combined_filepath_tsankov))
#   return(cbind(combined_count_overlaps, shendure_combined_count_overlaps,
#                tsankov_combined_count_overlaps))
# }

# colnames(combined_count_overlaps)[grep("Skin", colnames(combined_count_overlaps))]
# sum(combined_count_overlaps[grep("Skin T lymphocyte 2 \\(CD4\\+\\)", 
#                              colnames(combined_count_overlaps)), ])
# sum(combined_count_overlaps[grep("Skin Sun Exposed Macrophage \\(General,Alveolar\\)", 
#                                  colnames(combined_count_overlaps)), ])
# grep("Skin Fibroblast \\(Epithelial\\)", colnames(combined_count_overlaps))
# colnames(combined_count_overlaps)[grep("Skin", 
#                                        colnames(combined_count_overlaps))]

# cell_types = c("Skin Sun Exposed Melanocyte", 
#                "Skin Melanocyte",
#                "Skin Sun Exposed Fibroblast (Epithelial)",
#                "Skin Fibroblast (Epithelial)",
#                "Skin Keratinocyte 1", 
#                "Skin Sun Exposed Keratinocyte 1",
#                "Skin T Lymphocyte 1 (CD8+)", 
#                "Skin Sun Exposed T Lymphocyte 1 (CD8+)",
#                "Skin T lymphocyte 2 (CD4+)",
#                "Skin Sun Exposed T lymphocyte 2 (CD4+)",
#                "Skin Macrophage (General,Alveolar)",
#                "Skin Sun Exposed Macrophage (General,Alveolar)")
# 
# prep_boxplots_per_cancer_type(combined_count_overlaps,
#                               mut_count_data,
#                               "Skin.Melanoma",
#                               cell_types,
#                               1000,
#                               T,
#                               "skin_log_num_frags_vs_correlation.png")

# prep_boxplots(combined_count_overlaps, "Skin.Melanoma", cell_types, 200000, F,
#               "num_frags_vs_correlation.png")

# combined_counts_overlaps_all_scATAC_data = combine_scATAC(combined_filepath,
#                                                           combined_filepath_shendure,
#                                                           combined_filepath_tsankov)
combined_count_overlaps = t(readRDS(combined_filepath))
shendure_combined_count_overlaps = t(readRDS(combined_filepath_shendure))
tsankov_combined_count_overlaps = t(readRDS(combined_filepath_tsankov))
combined_counts_overlaps_all_scATAC_data = cbind(combined_count_overlaps, 
                                                 shendure_combined_count_overlaps,
                                                 tsankov_combined_count_overlaps)

# # RECALL that Shendure and Tsankov do NOT have Melanocytes
# colnames(combined_count_overlaps)[grep("Lung", 
#                                        colnames(combined_count_overlaps))]
# colnames(shendure_combined_count_overlaps)[grep("Lung", 
#                                     colnames(shendure_combined_count_overlaps))]
# colnames(tsankov_combined_count_overlaps)[grep("Lung", 
#                                      colnames(tsankov_combined_count_overlaps))]

lung_cell_types = c("Lung Alveolar Type 2 (AT2) Cell",
                    "Lung Bronchiolar and alveolar epithelial cells",
                    "Distal Lung AT2")

# sum(combined_counts_overlaps_all_scATAC_data[grep("Lung Bronchiolar and alveolar epithelial cells",
#                          colnames(combined_counts_overlaps_all_scATAC_data)), ])
lung_tsse_filtered_count_overlaps = c()
tsse_filtered_correlations = c()
for (file in mixedsort(list.files("processed_data/count_overlap_data/tsse_filtered/lung",
                        full.names = TRUE))) {
  count_overlaps = readRDS(file)
  tsse_filtered_correlations = append(tsse_filtered_correlations,
                                      cor(count_overlaps, 
                                      mut_count_data[, "Lung.AdenoCA"], 
                                      use="complete"))
}



prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
                              mut_count_data,
                              "Lung.AdenoCA",
                              lung_cell_types,
                              1000,
                              T,
                              "test_lung_log_num_frags_vs_correlation.png",
                              tsse_filtered_correlations)

skin_cell_types = c("Skin Sun Exposed Melanocyte", 
                    "Skin Fibroblast (Epithelial)")
skin_tsse_filtered_count_overlaps = c()
tsse_filtered_correlations = c()
for (file in mixedsort(list.files("processed_data/count_overlap_data/tsse_filtered/skin",
                                  full.names = TRUE))) {
  count_overlaps = readRDS(file)
  tsse_filtered_correlations = append(tsse_filtered_correlations,
                                      cor(count_overlaps, 
                                          mut_count_data[, "Skin.Melanoma"], 
                                          use="complete"))
}

prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
                              mut_count_data,
                              "Skin.Melanoma",
                              skin_cell_types,
                              1000,
                              T,
                              "test_skin_log_num_frags_vs_correlation.png",
                              tsse_filtered_correlations)
# # TSS
# get_filtered_metadata <- function(metadata, life_stage, tis, cell_types) {
#   cell_types = paste(cell_types, collapse="|")
#   filtered_metadata = metadata %>% 
#                       filter(grepl(life_stage, Life.stage)) %>%
#                       filter(grepl(tis, tissue, ignore.case=T)) %>%
#                       filter(grepl(cell_types, 
#                                    cell.type, 
#                                    ignore.case=T)) 
#   return(filtered_metadata)
# }
# 
# add_fragment_counts_to_metadata <- function(metadata, combined_count_overlaps) {
#   where_to_substr = unlist(lapply("_SM", regexpr, metadata[["tissue"]]))
#   tissue_name = unname(unlist(
#            lapply(metadata["tissue"], substr, 1, where_to_substr - 1)))
#   tissue_name = str_to_title(gsub("\\_", " ", tissue_name))
#   metadata = metadata %>% mutate(full_cell_type_name = paste(tissue_name, 
#                                                              cell.type))
#   count_overlaps_idx_per_cell_type = match(metadata[["full_cell_type_name"]], 
#                                            colnames(combined_count_overlaps))
#   uniq_count_overlap_idxs = unique(count_overlaps_idx_per_cell_type)
#   uniq_counts = apply(combined_count_overlaps[, uniq_count_overlap_idxs], 2, 
#                       sum)
#   names(uniq_counts) = uniq_count_overlap_idxs
#   count_overlap_idx_per_cell_type = match(count_overlaps_idx_per_cell_type, 
#                                           names(uniq_counts))
#   count_overlaps = unname(uniq_counts[count_overlap_idx_per_cell_type])
#   metadata["total_count_overlaps"] = count_overlaps
#   return(metadata)
# }
# 
# plot_num_frag_vs_tss_boxplot <- function(metadata, save_filename) {
#   ggplot(metadata) +
#     geom_boxplot(aes(x=factor(total_count_overlaps, levels=unique(total_count_overlaps)), 
#                      y=tsse, color=full_cell_type_name)) +
#     xlab("Num Fragments") +
#     theme(axis.text.x=element_text(size=10, angle = 90)) +
#     scale_x_discrete(limits = as.factor(sort(as.integer(
#       unique(metadata["total_count_overlaps"]) %>% pull))))
#   ggsave(paste("figures", save_filename, sep="/"), width = 20, height = 12)
# }
# 
# metadata = as_tibble(read.table("raw_dir/metadata/GSE184462_metadata.tsv", 
#                                 sep="\t",
#                       header=T))
# 
# filtered_metadata = get_filtered_metadata(metadata, "Adult", "skin", 
#                                           c("melanocyte",
#                                             "keratinocyte",
#                                             "fibroblast"))
# filtered_metadata = add_fragment_counts_to_metadata(filtered_metadata, 
#                                                     combined_count_overlaps)
# 
# plot_num_frag_vs_tss_boxplot(filtered_metadata, 
#                              "skin_num_frag_vs_tss.png")


# metadata = readRDS("processed_data/count_overlap_data/combined_count_overlaps/count_filter_1_combined_count_overlaps_metadata.rds")


# Num Cells
# melanocyte_n_cells = metadata[metadata["tissue_name"] == 
#                              "Skin" & metadata["cell_type"] == 
#                              "Melanocyte", "num_cells"]
