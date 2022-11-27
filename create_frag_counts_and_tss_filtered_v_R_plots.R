library(optparse)
library(gtools)
library(parallel)
library(tidyverse)
library(data.table)
source("utils.R")
options(scipen=999)

option_list <- list( 
  make_option("--cancer_type", type="character"),
  make_option("--boxplot_cell_types", type="character"),
  make_option("--tissue_for_tsse_filtered_cell_types", type="character",
              help="format: comma separated [Dataset]-[Tissue]"),
  make_option("--tsse_filtered_cell_types", type="character", 
              help="format: comma separated [CellType1]-[CellType2]-etc..., 
              order must correspond to tissue_for_tsse_filtered_cell_types"),
  make_option("--plot_filename", type="character"),
  make_option("--plot_x_ticks", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args =
#                     c("--cancer_type=Eso.AdenoCA",
#                       "--boxplot_cell_types=Stomach Goblet cells (SH)-Stomach Foveolar Cell (BR)-Colon Transverse Colon Epithelial Cell 2 (BR)-Esophagus Muscularis Foveolar Cell (BR)-Small Intestine Small Intestinal Enterocyte (BR)",
#                       "--tissue_for_tsse_filtered_cell_types=Shendure-Stomach,Bing Ren-Stomach,Bing Ren-Colon Transverse,Bing Ren-Small Intestine",
#                       "--tsse_filtered_cell_types=Goblet cells,Foveolar Cell,Colon Epithelial Cell 2,Small Intestinal Enterocyte",
#                       "--plot_filename=esophagus_adenoca_num_frags_vs_correlation.png",
#                       "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,600000"))

# Melanoma
# args = parse_args(OptionParser(option_list=option_list), args =
#                     c("--cancer_type=Skin.Melanoma",
#                       "--boxplot_cell_types=Skin Sun Exposed Melanocyte (BR)-Heart Schwann cells (SH)-Muscle Type II Skeletal Myocyte (BR)-Stomach Stromal cells (SH)-Lung Ciliated epithelial cells (SH)",
#                       "--tissue_for_tsse_filtered_cell_types=Bing Ren-Skin,Shendure-Heart,Bing Ren-Muscle,Shendure-Stomach,Shendure-Lung",
#                       "--tsse_filtered_cell_types=Melanocyte,Schwann cells,Type II Skeletal Myocyte,Stromal cells,Ciliated epithelial cells",
#                       "--plot_filename=melanoma_num_frags_vs_correlation.png",
#                       "--plot_x_ticks=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"))

#args = parse_args(OptionParser(option_list=option_list), args =
#                    c("--cancer_type=Lung.AdenoCA",
#                      "--boxplot_cell_types=Distal Lung AT2 (TS)-Lung Alveolar Type 2 (AT2) Cell (BR)-Lung Bronchiolar and alveolar epithelial cells (SH)",
#                      "--tissue_for_tsse_filtered_cell_types=Tsankov-Lung,Bing Ren-Lung,Shendure-Lung",
#                      "--tsse_filtered_cell_types=AT2,Alveolar Type 2 (AT2) Cell,Bronchiolar and alveolar epithelial cells",
#                      "--plot_filename=lung_adenoca_at2_only_frags_vs_correlation.png",
#                      "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,8000000,10000000,15000000,20000000"))

# args = parse_args(OptionParser(option_list=option_list), args =
#                    c("--cancer_type=Eso.AdenoCA",
#                      "--boxplot_cell_types=Stomach Goblet cells (SH)-Stomach Foveolar Cell (BR)-Colon Transverse Colon Epithelial Cell 2 (BR)-Esophagus Muscularis Foveolar Cell (BR)-Small Intestine Small Intestinal Enterocyte (BR)",
#                      "--tissue_for_tsse_filtered_cell_types=Shendure-Stomach,Bing Ren-Stomach,Bing Ren-Colon Transverse,Bing Ren-Small Intestine",
#                      "--tsse_filtered_cell_types=Goblet cells,Foveolar Cell,Colon Epithelial Cell 2,Small Intestinal Enterocyte",
#                      "--plot_filename=esophagus_adenoca_num_frags_vs_correlation.png",
#                      "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000"))

# args = parse_args(OptionParser(option_list=option_list), args =
#                    c("--cancer_type=Breast.AdenoCA",
#                      "--boxplot_cell_types=Mammary Tissue Basal Epithelial (Mammary) (BR)-Stomach Goblet cells (SH)-Mammary Tissue Mammary Luminal Epithelial Cell 1 (BR)-Mammary Tissue Mammary Luminal Epithelial Cell 2 (BR)-Esophagus Mucosa Airway Goblet Cell (BR)",
#                      "--tissue_for_tsse_filtered_cell_types=Bing Ren-Mammary Tissue,Shendure-Stomach,Bing Ren-Mammary Tissue,Bing Ren-Mammary Tissue,Bing Ren-Esophagus Mucosa",
#                      "--tsse_filtered_cell_types=Basal Epithelial (Mammary),Goblet cells,Mammary Luminal Epithelial Cell 1,Mammary Luminal Epithelial Cell 2,Airway Goblet Cell",
#                      "--plot_filename=breast_adenoca_num_frags_vs_correlation.png",
#                      "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000,15000000,20000000"))

# args = parse_args(OptionParser(option_list=option_list), args =
#                    c("--cancer_type=ColoRect.AdenoCA",
#                      "--boxplot_cell_types=Intestine Intestinal epithelial cells (SH)-Stomach Goblet cells (SH)-Placenta PAEP_MECOM positive cells (SH)-Colon Transverse Colon Epithelial Cell 2 (BR)-Liver Erythroblasts (SH)",
#                      "--tissue_for_tsse_filtered_cell_types=Shendure-Intestine,Shendure-Stomach,Shendure-Placenta,Bing Ren-Colon Transverse,Bing Ren-Liver",
#                      "--tsse_filtered_cell_types=Intestinal epithelial cells,Goblet cells,PAEP_MECOM positive cells,Colon Epithelial Cell 2,Erythroblasts",
#                      "--plot_filename=colorect_adenoca_num_frags_vs_correlation.png",
#                      "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000"))

# args = parse_args(OptionParser(option_list=option_list), args =
#                    c("--cancer_type=CNS.GBM",
#                      "--boxplot_cell_types=nerve_tibial Schwann Cell (General) (BR)/heart_lv Cardiac Pericyte 1 (BR)/muscle Type I Skeletal Myocyte (BR)/lung Club Cell (BR)/heart_lv Cardiac Pericyte 2 (BR)/heart_atrial_appendage Cardiac Pericyte 1 (BR)/heart_lv Endothelial Cell (General) 1 (BR)/lung Macrophage (General) (BR)/adipose_omentum Macrophage (General,Alveolar) (BR)/colon_transverse Smooth Muscle (General) (BR)",
#                      "--tissue_for_tsse_filtered_cell_types=Bing Ren-Nerve Tibial,Bing Ren-Heart Lv,Bing Ren-Muscle,Bing Ren-Lung,Bing Ren-Heart Lv,Bing Ren-Heart Atrial Appendage,Bing Ren-Heart Lv,Bing Ren-Lung,Bing Ren-Adipose Omentum,Bing Ren-Colon Transverse",
#                      "--tsse_filtered_cell_types=Schwann Cell (General)/Cardiac Pericyte 1/Type I Skeletal Myocyte/Club Cell/Cardiac Pericyte 2/Cardiac Pericyte 1/Endothelial Cell (General) 1/Macrophage (General)/Macrophage (General,Alveolar)/Smooth Muscle (General)",
#                      "--plot_filename=cns_gbm_BR_only_with_top_20_bottom_5_feats_num_frags_vs_correlation.png",
#                      "--plot_x_tick=1000,10000,50000,100000,150000,250000,300000,400000,500000,1000000,1500000,2000000,5000000,10000000"))

cancer_type = args$cancer_type
boxplot_cell_types = unlist(strsplit(args$boxplot_cell_types, split = "/"))
tissue_dataset_for_tsse_filtered_cell_types = unlist(strsplit(args$tissue_for_tsse_filtered_cell_types, 
                                              split = ","))
tsse_filtered_cell_types = unlist(strsplit(args$tsse_filtered_cell_types, 
                                           split = "/"))
plot_filename = args$plot_filename
plot_x_ticks = as.integer(unlist(strsplit(args$plot_x_ticks, split = ",")))

CELL_NUMBER_FILTER = 1
load('/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

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

# create_plot_x_ticks <- function(lower, upper, log=T) {
#   if (log) {
#     lower = log(lower, base=2)
#     upper = log(upper, base=2)
#   }
#   num_fragments_to_subsample <- seq(lower, upper, (upper - lower) / 9)
#   # num_fragments_to_subsample[2:(length(num_fragments_to_subsample) - 1)] = 
#   #   num_fragments_to_subsample[2:(length(num_fragments_to_subsample) 
#   #                                            - 1)]
#   if (!log) {
#     num_fragments_to_subsample[2:(length(num_fragments_to_subsample) 
#                                   - 1)] =
#       as.integer(num_fragments_to_subsample[2:(length(num_fragments_to_subsample)
#                                                - 1)])
#   }
#   return(num_fragments_to_subsample)
# }

run_simulations_per_num_fragments <- function(num_fragments) {
  # if (log) {
  #   num_fragments = as.integer(2^num_fragments)
  # }
  samples_10 = mclapply(rep(i, 10), run_one_simulation, 
                        cell_type_count_overlaps,
                        sampling_vec, mutations,
                        mc.cores=8)
  return(unlist(samples_10))
}

run_subsampling_simulations <- function(num_fragments_to_subsample, 
                                        cell_type_count_overlaps, 
                                        sampling_vec, mutations, cell_type) {
  subsampling_simulations <- list()
  list_idx = 1
  set.seed(42)
  
  # subsampling_simulations = mclapply(num_fragments_to_subsample, 
  #                                    run_simulations_per_num_fragments,
  #                                    log, mc.cores=8)
  for (i in num_fragments_to_subsample) {
    ptm <- proc.time()
    # if (log) {
    #   i = as.integer(2^i)
    # }
    print(paste("subsampling", i, 
                "fragments for", paste0(cell_type), sep=" "))
    samples_10 = mclapply(rep(i, 10), run_one_simulation, 
                          cell_type_count_overlaps,
                          sampling_vec, mutations,
                          mc.cores=8)
    subsampling_simulations[[list_idx]] = unlist(samples_10)
    list_idx = list_idx + 1
    print(paste("time to subsample", i, 
                "for", paste0(cell_type, ":"), round((proc.time() - ptm)[["elapsed"]]),
                "seconds"), sep=" ")
  }
  names(subsampling_simulations) = num_fragments_to_subsample
  
  # if (log) {
  #   names(subsampling_simulations) = as.integer(unlist(lapply(2, `^`, 
  #                                               num_fragments_to_subsample)))
  # }
  subsampling_simulations = as_tibble(subsampling_simulations)
  return(subsampling_simulations)
}

add_escape_helper <- function(char_to_escape, cell_type) {
  escaped = paste0("\\", char_to_escape)
  left_idx = unlist(gregexpr(escaped, cell_type))
  add = 0:(length(left_idx) - 1)
  left_idx = left_idx + add
  cell_type_to_grep = cell_type
  for (idx in left_idx) {
    left_substr = substr(cell_type_to_grep, 1, idx - 1)
    right_substr = substr(cell_type_to_grep, idx, nchar(cell_type_to_grep))
    cell_type_to_grep = paste0(left_substr, "\\", right_substr)
  }
  return(cell_type_to_grep)
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
  if (grepl("\\?", cell_type)) {
    cell_type = add_escape_helper("?", cell_type)
  }
  return(cell_type)
}

plot_and_save_boxplots <- function(correlations_long, cell_types, 
                                   plot_filename, 
                                   correlations_for_tsse_filtered_cells=NULL) {
  switch = F
  if (!is.null(correlations_for_tsse_filtered_cells)) {
    # x_for_filtered_correlations = unique(correlations_long %>% 
    #                                        pull("num_fragments"))
    # correlations_tsse_filtered = 
    #   as.list(as.data.frame(t(correlations_for_tsse_filtered_cells)))
    
    # correlations_tsse_filtered = lapply(correlations_tsse_filtered,
    #                                     function(x) x[!is.na(x)])
    # x_for_filtered_correlations = lapply(correlations_tsse_filtered,
    #                                      function(x) 
    #                                    x_for_filtered_correlations[1:length(x)])
    
    # levels_for_filtered_correlations = x_for_filtered_correlations
    cell_type_for_tsse_filtered = lapply(correlations_for_tsse_filtered_cells, 
                                         rownames)
    dfs = lapply(correlations_for_tsse_filtered_cells, 
                 as_tibble)
    df = rbindlist(dfs)
    df[, "cell_type"] = unlist(cell_type_for_tsse_filtered)
    df = pivot_longer(df, -cell_type)
    colnames(df)[2] = "frag_counts"
    df["frag_counts"] = as.integer(df[["frag_counts"]])
    # lengths = lapply(x_for_filtered_correlations, length)
    
    # repeat_cell_types <- function(i, cell_type_for_tsse_filtered, lengths) {
    #   return(rep(cell_type_for_tsse_filtered[[i]], each=lengths[[i]]))
    # }
    # 
    # cell_type_for_tsse_filtered = lapply(seq_along(cell_type_for_tsse_filtered),
    #                                      repeat_cell_types,
    #                                      cell_type_for_tsse_filtered,
    #                                      lengths)
    # 
    # create_cell_type_and_frag_counts_columns <- function(i, df, 
    #                                                      cell_type_for_tsse_filtered,
    #                                                      x_for_filtered_correlations) {
    #   df[[i]]["frag_counts"] = x_for_filtered_correlations[[i]]
    #   df[[i]]["cell_type"] = cell_type_for_tsse_filtered[[i]]
    #   
    #   return(df[[i]])
    # }
    # 
    # dfs = lapply(seq_along(dfs), create_cell_type_and_frag_counts_columns, 
    #              dfs, cell_type_for_tsse_filtered, x_for_filtered_correlations)
    # df = rbindlist(dfs)
    switch = T
  }
  ggplot() +
    {if (switch) geom_point(data=df, aes(x=factor(frag_counts, 
                                                  levels=unique(frag_counts)),
                                         y=value,
                                         fill=cell_type), shape=21,
                            size=3)} +
    geom_boxplot(data=correlations_long,
                 aes(x=factor(num_fragments, levels=unique(num_fragments)), 
                     y=correlation, color=cell_type)) +
    theme(axis.text.x=element_text(size=10, angle=90)) +
    xlab("Num Fragments") +
    ylab("Correlation") +
    labs(color="Cell type", fill = "Cell type (TSS filtered)") +
    scale_x_discrete(limits = as.factor(sort(as.integer(
      unique(correlations_long["num_fragments"]) %>% pull))))#+ 
  ggsave(paste("figures", plot_filename, sep="/"), width = 20, height = 12)
}

remove_bigger_than_curr_fragments <- function(plot_x_ticks,
                                              curr_num_fragments) {
  plot_x_ticks_reduced <- c()
  for (num_fragments in plot_x_ticks) {
    # if (log) {
    #   if (as.integer(2^num_fragments) <= curr_num_fragments) {
    #     plot_x_ticks_reduced <- append(plot_x_ticks_reduced,
    #                                                   num_fragments)
    #   }
    # }
    if (num_fragments <= curr_num_fragments) {
      plot_x_ticks_reduced <- append(plot_x_ticks_reduced,
                                                    num_fragments)
    }
  }
  return(plot_x_ticks_reduced)
}

get_long_correlations_per_cell_type <- function(i,
                                                cell_types, 
                                                cell_types_total_fragments,
                                                plot_x_ticks, 
                                                mutations, 
                                                combined_count_overlaps) {
  cell_type_total_fragments = cell_types_total_fragments[[i]]
  cell_type_count_overlaps = combined_count_overlaps[, i]
  cell_type_plot_x_ticks = remove_bigger_than_curr_fragments(
    plot_x_ticks[[i]],
    cell_type_total_fragments)
  if (is.null(cell_type_plot_x_ticks)) {
    return(tibble())
  }
  # upper = cell_type_total_fragments
  sampling_vec = rep(names(cell_type_count_overlaps), 
                     cell_type_count_overlaps)
  # plot_x_ticks <- create_plot_x_ticks(lower, 
  #                                                                   upper, 
  #                                                                   log)
  subsampling_simulations <- run_subsampling_simulations(cell_type_plot_x_ticks, 
                                                         cell_type_count_overlaps, 
                                                         sampling_vec, 
                                                         mutations,
                                                         cell_types[[i]])
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

get_missing_fragment_counts_per_cell_type <- function(type, correlations_long,
                                                      plot_x_ticks) {
  existing_frag_counts = as.integer(unique(correlations_long %>% 
                                  filter(cell_type == type) %>%
                                  pull(num_fragments)))
  missing_fragment_counts = plot_x_ticks[!(plot_x_ticks %in%
                                           existing_frag_counts)]
  return(missing_fragment_counts)
}

get_correlations_from_missing_fragment_counts <- function(correlations_long,
                                                          plot_x_ticks,
                                                          cell_types,
                                                          cell_type_count_overlaps,
                                                          total_fragments_per_cell_type,
                                                          mutations) {
  missing_fragment_counts_per_cell_type = mclapply(cell_types, 
                                                   get_missing_fragment_counts_per_cell_type,
                                                   correlations_long,
                                                   plot_x_ticks,
                                                   mc.cores=8)
  correlations_long = mclapply(seq_along(cell_types),
                               get_long_correlations_per_cell_type,
                               cell_types,
                               total_fragments_per_cell_type,
                               missing_fragment_counts_per_cell_type,
                               mutations,
                               cell_type_count_overlaps,
                               mc.cores=8)
  return(correlations_long)
}

add_na_correlations <- function(correlations, plot_x_ticks) {
  for (frag_counts in plot_x_ticks) {
    if (!(frag_counts %in% colnames(correlations))) {
      correlations = correlations %>% add_column(frag_counts = NA)
      colnames(correlations)[length(colnames(correlations))] = frag_counts
    }
  }
  return(correlations)
}

prep_boxplots_per_cancer_type <- function(combined_count_overlaps, 
                                          mut_count_data, 
                                          cancer_type, 
                                          cell_types, 
                                          plot_x_ticks,
                                          plot_filename, 
                                          correlations_for_tsse_filtered_cells=NULL) {
  cell_type_for_grep = lapply(cell_types, add_escape_if_necessary)
  cell_type_col_idx = lapply(cell_type_for_grep, grep,
                             colnames(combined_count_overlaps), ignore.case=T)
  cell_types_count_overlaps = combined_count_overlaps[, 
                                                      unlist(cell_type_col_idx)]
  cell_types_total_fragments = apply(cell_types_count_overlaps, 2, sum)
  mutations <- mut_count_data[cancer_type]
  correlations_filename = paste0(cancer_type, ".rds")
  correlations_filepath = paste("processed_data", "subsampled_correlations",
                                correlations_filename, sep="/")
  
  if (!file.exists(correlations_filepath)) {
    # plot_x_ticks <- create_plot_x_ticks(lower, 
    #                                             max(cell_types_total_fragments), 
    #                                             log)
    
    correlations_long = mclapply(seq_along(cell_types),
                                 get_long_correlations_per_cell_type,
                                 colnames(cell_types_count_overlaps),
                                 cell_types_total_fragments,
                                 rep(list(plot_x_ticks), 
                                     length(cell_types)),
                                 mutations,
                                 cell_types_count_overlaps,
                                 mc.cores=8)
    correlations_long = bind_rows(correlations_long)
  }
  else {
    correlations_long = readRDS(correlations_filepath)
    correlations_long_missing = get_correlations_from_missing_fragment_counts(
      correlations_long, 
      plot_x_ticks,
      cell_types,
      cell_types_count_overlaps,
      cell_types_total_fragments,
      mutations)
    correlations_long_missing = bind_rows(correlations_long_missing)
    correlations_long = rbind(correlations_long, correlations_long_missing)
  }
  saveRDS(correlations_long, correlations_filepath)
  correlations_long["num_fragments"] = as.integer(correlations_long[["num_fragments"]])
  correlations_long = correlations_long %>% 
                      filter(num_fragments %in% plot_x_ticks) %>%
                      filter(cell_type %in% cell_types)
  correlations_for_tsse_filtered_cells = lapply(correlations_for_tsse_filtered_cells,
                                                add_na_correlations, 
                                                plot_x_ticks)
  plot_and_save_boxplots(correlations_long, colnames(cell_types_count_overlaps),
                         plot_filename, correlations_for_tsse_filtered_cells)
}

combine_tsse_filtered_count_overlaps_into_correlation_df <- function(folder_path,
                                                                     cancer_type,
                                                                     plot_x_ticks,
                                                                     current_tsse_cell_types) {
  cell_type_for_grep = lapply(current_tsse_cell_types, add_escape_if_necessary)
  count = 1
  tsse_filtered_correlations = c()
  for (file in mixedsort(list.files(folder_path,
                                    full.names = TRUE))) {
    num_frags = as.integer(unlist(strsplit(tail(unlist(strsplit(file, "_")), 1), 
                                "[.]"))[1])
    if (num_frags %in% plot_x_ticks) {
      num_fragments = tail(unlist(strsplit(file, split="/")), 1)
      num_fragments = tail(unlist(strsplit(num_fragments, split="_")), 1)
      num_fragments = head(unlist(strsplit(num_fragments, split="[.]")), 1)
      
      count_overlaps = readRDS(file)
      cell_type_col_idx = lapply(cell_type_for_grep, grep,
                                 colnames(count_overlaps), ignore.case=T)

      count_overlaps = count_overlaps[, unlist(cell_type_col_idx), drop=FALSE]
      if (length(count_overlaps) == 0) {
        corrs = data.frame()
      }
      else {
      corrs = data.frame(cor(count_overlaps,
                             mut_count_data[, cancer_type],
                             use="complete"))
      }
      for (cell_type in rownames(tsse_filtered_correlations)) {
        if (!(cell_type %in% rownames(corrs))) {
          corrs[cell_type, 1] = c(NA)
        }
      }
      colnames(corrs) = "correlation"
      if (count == 1) {
        tsse_filtered_correlations = corrs
      }
      else {
        tsse_filtered_correlations = cbind(tsse_filtered_correlations, corrs)
      }
      colnames(tsse_filtered_correlations)[ncol(tsse_filtered_correlations)] = 
        num_fragments
      count = count + 1
    }
  }
  return(tsse_filtered_correlations)
}

# combine_scATAC <- function(combined_filepath, combined_filepath_shendure,
#                            combined_filepath_tsankov) {
#   combined_count_overlaps = t(readRDS(combined_filepath))
#   shendure_combined_count_overlaps = t(readRDS(combined_filepath_shendure))
#   tsankov_combined_count_overlaps = t(readRDS(combined_filepath_tsankov))
#   return(cbind(combined_count_overlaps, shendure_combined_count_overlaps,
#                tsankov_combined_count_overlaps))
# }


# RECALL that Shendure and Tsankov do NOT have Melanocytes
# colnames(combined_count_overlaps)[grep("Lung",
#                                         colnames(combined_count_overlaps))]
# colnames(shendure_combined_count_overlaps)[grep("Lung",
#                                      colnames(shendure_combined_count_overlaps))]
# colnames(tsankov_combined_count_overlaps)[grep("Lung",
#                                       colnames(tsankov_combined_count_overlaps))]

# lung_cell_types = c("Lung Alveolar Type 2 (AT2) Cell",
#                      "Lung Bronchiolar and alveolar epithelial cells",
#                      "Distal Lung AT2")
# 
# # lung_tsse_filtered_count_overlaps = c()
# # tsse_filtered_correlations = c()
# # for (file in mixedsort(list.files("processed_data/count_overlap_data/tsse_filtered/lung",
# #                          full.names = TRUE))) {
# #    count_overlaps = readRDS(file)
# #    tsse_filtered_correlations = append(tsse_filtered_correlations,
# #                                        cor(count_overlaps,
# #                                        mut_count_data[, "Lung.AdenoCA"],
# #                                        use="complete"))
# # }
# lung_tsse_filtered_correlations =
#   combine_tsse_filtered_count_overlaps_into_correlation_df(
#     "processed_data/count_overlap_data/tsse_filtered/lung",
#     "Lung.AdenoCA")
# 
# rownames(lung_tsse_filtered_correlations) =
#   paste("Lung", rownames(lung_tsse_filtered_correlations))
# 
# 
# tsse_filtered_correlations = list(lung_tsse_filtered_correlations)
# prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
#                                mut_count_data,
#                                "Lung.AdenoCA",
#                                lung_cell_types,
#                                1000,
#                                T,
#                                "lung_adenoca_log_num_frags_vs_correlation.png",
#                                tsse_filtered_correlations)

mut = load_mutation_data()
mut_count_data = get_mutation_df_all_cancers()
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

combined_count_overlaps = t(readRDS(combined_filepath))
colnames(combined_count_overlaps) = paste(colnames(combined_count_overlaps), 
                                          "(BR)")
shendure_combined_count_overlaps = t(readRDS(combined_filepath_shendure))
colnames(shendure_combined_count_overlaps) = paste(colnames(shendure_combined_count_overlaps), 
                                          "(SH)")
tsankov_combined_count_overlaps = t(readRDS(combined_filepath_tsankov))
colnames(tsankov_combined_count_overlaps) = paste(colnames(tsankov_combined_count_overlaps), 
                                                  "(TS)")
combined_counts_overlaps_all_scATAC_data = cbind(combined_count_overlaps, 
                                                 shendure_combined_count_overlaps,
                                                 tsankov_combined_count_overlaps)
tsse_filtered_correlations = list()
idx = 1
for (dataset_tissue in tissue_dataset_for_tsse_filtered_cell_types) {
  dataset = unlist(strsplit(dataset_tissue, split = "-"))[1]
  tissue = unlist(strsplit(dataset_tissue, split = "-"))[2]
  dataset_dir = tolower(dataset)
  dataset_dir = gsub(" ", "_", dataset_dir)
  tissue_dir = tolower(tissue)
  tissue_dir = gsub(" ", "_", tissue_dir)
  path = paste("processed_data/count_overlap_data/tsse_filtered",
               dataset_dir, tissue_dir, sep="/")
  current_tsse_cell_types = tsse_filtered_cell_types[idx]
  current_tsse_cell_types = unlist(strsplit(current_tsse_cell_types, 
                                            split = ";"))
  correlations =
    combine_tsse_filtered_count_overlaps_into_correlation_df(
      path,
      cancer_type, 
      plot_x_ticks,
      current_tsse_cell_types)
  if (dataset == "Bing Ren") {
    dataset_extension = "(BR)"
  }
  else if (dataset == "Shendure") {
    dataset_extension = "(SH)"
  }
  else {
    dataset_extension = "(TS)"
  }
  
  rownames(correlations) =
    paste(tissue, rownames(correlations), dataset_extension)
  
  tsse_filtered_correlations = append(tsse_filtered_correlations,
                                      list(correlations))
  idx = idx + 1
}

prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
                              mut_count_data,
                              cancer_type,
                              boxplot_cell_types,
                              plot_x_ticks,
                              plot_filename,
                              tsse_filtered_correlations)


# #### ColoRect.AdenoCA ####
# cell_types = c("Colon Transverse Colon Epithelial Cell 2",
#                 "Colon Transverse T Lymphocyte 1 (CD8+)",
#                 "Mammary Tissue Basal Epithelial (Mammary)",
#                 "Mammary Tissue Mammary Luminal Epithelial Cell 1",
#                 "Colon Transverse Colonic Goblet Cell")
# 
# colon_transverse_tsse_filtered_correlations =
#   combine_tsse_filtered_count_overlaps_into_correlation_df(
#     "processed_data/count_overlap_data/tsse_filtered/colon_transverse",
#     "ColoRect.AdenoCA")
# 
# rownames(colon_transverse_tsse_filtered_correlations) =
#   paste("Colon Transverse",
#         rownames(colon_transverse_tsse_filtered_correlations))
# 
# mammary_tissue_tsse_filtered_correlations =
#   combine_tsse_filtered_count_overlaps_into_correlation_df(
#     "processed_data/count_overlap_data/tsse_filtered/mammary_tissue",
#     "ColoRect.AdenoCA")
# 
# rownames(mammary_tissue_tsse_filtered_correlations) =
#   paste("Mammary Tissue", rownames(mammary_tissue_tsse_filtered_correlations))
# 
# tsse_filtered_correlations = list(colon_transverse_tsse_filtered_correlations,
#                                   mammary_tissue_tsse_filtered_correlations)
# 
# prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
#                               mut_count_data,
#                               "ColoRect.AdenoCA",
#                               cell_types,
#                               1000,
#                               T,
#                               "colorect_adenoca_log_num_frags_vs_correlation.png",
#                               tsse_filtered_correlations)

# #### Breast AdenoCA ####
# cell_types = c("Mammary Tissue Mammary Luminal Epithelial Cell 1",
#                "Mammary Tissue Basal Epithelial (Mammary)",
#                "Esophagus Mucosa Airway Goblet Cell",
#                "Mammary Tissue Mammary Luminal Epithelial Cell 2",
#                "Heart Atrial Appendage Cardiac Pericyte 1")
# 
# mammary_tissue_tsse_filtered_correlations =
#   combine_tsse_filtered_count_overlaps_into_correlation_df(
#     "processed_data/count_overlap_data/tsse_filtered/mammary_tissue",
#     "Breast.AdenoCA")
# 
# prep_boxplots_per_cancer_type(combined_counts_overlaps_all_scATAC_data,
#                               mut_count_data,
#                               "Breast.AdenoCA",
#                               cell_types,
#                               1000,
#                               "breast_adenoca_log_num_frags_vs_correlation.png",
#                               mammary_tissue_tsse_filtered_correlations)
