library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tibble)
library(tidyr)
library(stringr)
library(exomeCopy)
source("utils.R")

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')
dir.create("figures") 

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

create_per_tissue_heatmap <- function(cor_df, n, i, frag_df) {
  tissue = n[i]
  tissue = gsub(" ", "_", tolower(tissue), )
  heatmap_filename = tissue
  correlations = cor_df[[i]]
  fragments_per_row = frag_df[[i]]
  
  path = paste("figures/per_tissue_heatmaps/", heatmap_filename, ".pdf", sep="")
  ha = rowAnnotation(fragment_counts = anno_barplot(fragments_per_row))
  save_heatmap_file(path, 40, 40, correlations, left_annotation=ha) 
  
  correlations_rowwise = t(scale(t(correlations)))
  correlations_columnwise = scale(correlations)
  row_heatmap_filename = paste(heatmap_filename, "row_normalized", sep="_")
  column_heatmap_filename = paste(heatmap_filename, "column_normalized", 
                                  sep="_")
  
  row_path = paste("figures/per_tissue_heatmaps/", row_heatmap_filename, ".pdf",
                   sep="")
  column_path = paste("figures/per_tissue_heatmaps/", column_heatmap_filename, 
                      ".pdf", sep="")
  
  save_heatmap_file(row_path, 40, 40, correlations_rowwise, left_annotation=ha)
  save_heatmap_file(column_path, 40, 40, correlations_columnwise, left_annotation=ha)
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

combine_count_overlaps_per_tissue <- function(f, tissue_specific_counts) {
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
  return(tissue_specific_counts)
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

#### Cell type Cell type correlation ####
# Expect cells of same type from different tissues to cluster together, rather
# than cells from same tissue to cluster 
cell_type_cell_type_corr = cor(combined_count_overlaps, combined_count_overlaps)
total_fragments = apply(combined_count_overlaps, 2, sum)
ha = rowAnnotation(fragment_counts = anno_barplot(log2(total_fragments)))
save_heatmap_file("figures/cell_type_cell_type_corr.pdf",
                  30, 30, cell_type_cell_type_corr, 
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6), left_annotation = ha)

save_heatmap_file("figures/cell_type_cell_type_corr_pearson.pdf",
                  30, 30, cell_type_cell_type_corr, 
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6),
                  left_annotation = ha, 
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson")

#### One heatmap, all tissues ####
correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
correlations = t(scale(t(correlations)))
correlations = scale(correlations)
save_heatmap_file("figures/all_cells_corr_double_normalized.pdf",
                  30, 30, correlations, row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6))

#### Tissue specific heatmaps ####
tissue_specific_counts_path = "processed_data/count_overlap_data/combined_count_overlaps/tissue_specific_counts_list.rds"
tissue_specific_counts = list()

if (!file.exists(tissue_specific_counts_path)) {
  pattern = "^count_filter_.*"
  for (f in list.files("processed_data/count_overlap_data/", pattern=pattern,
                          full.names = TRUE)) {
    tissue_specific_counts <- combine_count_overlaps_per_tissue(f, 
                                                         tissue_specific_counts)
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

#### Min correlation vs number of fragments for Fibroblast ####
cell_type = 'Fibroblast'
unsquashed_count_overlaps = readRDS("processed_data/count_overlap_data/combined_count_overlaps/count_filter_100_unsquashed_count_overlaps.rds")
unsquashed_tissue_names =  readRDS("processed_data/count_overlap_data/combined_count_overlaps/count_filter_100_unsquashed_tissue_names.rds")
total_fragments_per_sample = rowSums(unsquashed_count_overlaps)
unsquashed_correlations = cor(t(unsquashed_count_overlaps), mut_count_data, 
                              use = "complete.obs")
min_corr_per_sample = apply(unsquashed_correlations, 1, min)
df = data.frame(fragments=total_fragments_per_sample, 
                correlation=min_corr_per_sample,
                tissue=unsquashed_tissue_names)
df = df[grepl(cell_type, rownames(df), ignore.case = T), ]
corr = round(cor(df[,'fragments'], df[, 'correlation']), 2)
ggplot(df, aes(x=fragments, y=correlation, color=tissue)) +
  geom_point() +
  annotate("text", x=20000000, y=-0.2, label= toString(str_glue("R=", corr)))
  #stat_cor(method = "pearson",label.x = 5000000, label.y = -1.5)

#### Polak Count overlaps ####
polak_combined_count_overlaps = data.frame()
polak_combined_filepath = 
  "processed_data/polak_count_overlap_data/combined_count_overlaps/polak_combined_count_overlaps.rds"

if (!file.exists(polak_combined_filepath)) {
  polak_combined_count_overlaps = data.frame()
  for (f in list.files("processed_data/polak_count_overlap_data", 
                       full.names=TRUE, pattern=
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
  saveRDS(polak_combined_count_overlaps, polak_combined_filepath)
}

polak_combined_count_overlaps = t(readRDS(polak_combined_filepath))
correlations = cor(polak_combined_count_overlaps, mut_count_data, 
                   use="complete.obs")
correlations = scale(correlations)
save_heatmap_file("figures/polak_cell_type_to_cancer_correlations.pdf", 
                  10, 30, correlations, 
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6))

#### Our data vs Polak correlation ####
correlations = t(scale(t(cor(combined_count_overlaps, 
                             polak_combined_count_overlaps))))
correlations = scale(correlations)
save_heatmap_file("figures/our_data_vs_polak_corr.pdf", 10, 30, correlations,
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 6))

#### Cell number / Fragment counts vs R boxplots ####
correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
cell_types = rownames(cor(combined_count_overlaps, mut_count_data, use="complete.obs"))
correlations = as_tibble(correlations)
correlations = correlations %>%
               add_column(total_fragments=total_fragments) %>%
               add_column(cell_type=cell_types)

correlations_long = pivot_longer(correlations, cols = 
                                   -(total_fragments:cell_type),
                                 names_to = "cancer_type", 
                                 values_to = "correlation")
# correlations_long$total_fragments = 
  # factor(correlations_long$total_fragments, 
  #        levels = unique(correlations_long$total_fragments[order (correlations_long$total_fragments)]))
# geom_boxplot(aes(x=total_fragments, y=correlation)) +
# theme(axis.text.x=element_text(size=3, angle = 90)) +
# # xlab("Num Fragments") +
# scale_x_discrete(name="Num Fragments") +

medians <- correlations_long %>% 
  group_by(total_fragments) %>% 
  summarise(median(correlation))

colnames(medians)[2] <- "corr_medians"
correlations_long <- correlations_long %>% right_join(medians)

ggplot(correlations_long) +
  geom_boxplot(aes(x=as.factor(total_fragments), y=correlation)) +
  theme(axis.text.x=element_text(size=3, angle = 90)) +
  xlab("Num Fragments") 

correlations_long_no_outliers = correlations_long[correlations_long$corr_medians < -0.2, ]
medians = medians[medians$corr_medians < -0.2, ]
coefs <- coef(lm(corr_medians ~ total_fragments, data = medians))

ggplot(correlations_long_no_outliers) +
  geom_boxplot(aes(x=as.factor(total_fragments), y=correlation)) +
  theme(axis.text.x=element_text(size=3, angle = 90)) +
  xlab("Num Fragments") +
  geom_abline(intercept = coefs[1], slope = coefs[2])


  # stat_cor(aes(x=unique(total_fragments), y=unique(corr_medians)), method = "pearson")

# cor(medians$total_fragments, medians$corr_medians)

#+
  # geom_point(aes(x=as.factor(total_fragments), y=corr_medians)) +
  # geom_smooth(aes(x=total_fragments, y=corr_medians))

  # stat_cor(aes(x=unique(total_fragments), y=unique(corr_medians)), method = "pearson")
# 
# ggplot(medians, aes(x=total_fragments, y=corr_medians)) +
#   geom_point() +
#   stat_cor(method = "pearson")

# length(unique(correlations_long$total_fragments))
# length(unique(correlations_long$corr_medians))

# ggplot(correlations_long, aes(x=total_fragments, y=corr_medians)) +
#   stat_cor(method = "pearson")

# label.x = 5000000, label.y = 1.5)
#geom_smooth(aes(x=total_fragments, y=median))
# geom_boxplot(aes(x=total_fragments, y=correlation)) +
# theme(axis.text.x=element_text(size=3, angle = 90)) +
# # xlab("Num Fragments") +
# scale_x_discrete(name="Num Fragments") +