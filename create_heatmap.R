library(stringr)
library(ComplexHeatmap)

# proj_dir = '/ahg/regevdata/projects/ICA_Lung/Bruno/celloforigin_prj/'
# mut = readRDS(paste0(proj_dir, 'mutations.aggregated.PCAWG.RK.20181215.Rds'))
# load('/ahg/regevdata/projects/ICA_Lung/Bruno/celloforigin_prj/hg19.1Mb.ranges.Polak.Nature2015.RData')
# interval.ranges

mut = readRDS('mutations.aggregated.PCAWG.RK.20181215.Rds')
load('hg19.1Mb.ranges.Polak.Nature2015.RData')
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

create_per_tissue_heatmap <- function(cor_df, n, i, normalize) {
  tissue = n[i]
  tissue = gsub(" ", "_", tolower(tissue), )
  heatmap_filename = tissue
  correlations = cor_df[[i]]
  if (normalize) {
    correlations = t(scale(t(correlations)))
    heatmap_filename = paste(heatmap_filename, "normalized", sep="_")
  }

  path = paste("figures/per_tissue_heatmaps/", heatmap_filename, ".pdf", sep="")
  pdf(path)
  h = Heatmap(correlations, row_names_gp = gpar(fontsize = 8),
                            column_names_gp = gpar(fontsize = 6))
  print(h)
  dev.off()
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

pattern = paste("count_filter_", CELL_NUMBER_FILTER, "_", sep="")
if (!file.exists(combined_filepath)) {
  combined_count_overlaps = data.frame()
  for (file in list.files("count_overlap_data", pattern = pattern, 
                          full.names = TRUE)) {
    tissue_name <- get_tissue_name(file)
  	count_overlaps = readRDS(file)
  	count_overlaps = as.data.frame(do.call(rbind, count_overlaps), 
  	                               row.names = paste(tissue_name, 
  	                               names(count_overlaps)))
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
}

combined_filepath = paste("count_overlap_data/processed_count_overlaps/", 
                          "count_filter_",
                          100, 
                          "_combined_count_overlaps.rds", sep="")

combined_count_overlaps = t(readRDS(combined_filepath))

correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
correlations = t(scale(t(correlations)))
# correlations = scale(correlations)

write.csv(correlations, "all_corrs.csv")

#png(file="/home/mdanb/research/mount_sinai/bing_ren/corrs.png")
pdf ("/home/mdanb/research/mount_sinai/bing_ren/corrs.pdf", width = 10, height = 30)
h = Heatmap(correlations,
            row_names_gp = gpar(fontsize = 6),
            column_names_gp = gpar(fontsize = 6),height=18,
)
print(h)
dev.off()


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
per_tissue_cor = lapply(tissue_specific_counts, cor, mut_count_data, use="complete.obs")

dir.create("figures/per_tissue_heatmaps")
lapply(seq_along(per_tissue_cor), create_per_tissue_heatmap, 
       cor_df=per_tissue_cor, 
       n=names(per_tissue_cor),
       normalize=T)
lapply(seq_along(per_tissue_cor), create_per_tissue_heatmap, 
       cor_df=per_tissue_cor, 
       n=names(per_tissue_cor),
       normalize=F)

# Correlation analysis between num mutation and average correlation 
tot_mutations_per_cancer_type = colSums(mut_count_data, na.rm = T)
cor(t(correlations), tot_mutations_per_cancer_type)
mean(cor(t(correlations), tot_mutations_per_cancer_type))


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
