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

# tissue_specific_counts = list()
# for (file in list.files("count_overlap_data", pattern=pattern, 
#                         full.names = TRUE)) {
#   count_overlaps = readRDS(file)
#   tissue_name <- get_tissue_name(file)
# 
#   count_overlaps = as.data.frame(do.call(rbind, count_overlaps),
#                                  row.names = paste(tissue_name,
#                                                    names(count_overlaps)))
#   tissue_specific_counts[[tissue_name]] = append(
#                                           tissue_specific_counts[[tissue_name]],
#                                           count_overlaps)
# }


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


combined_count_overlaps = t(readRDS(combined_filepath))

correlations = cor(combined_count_overlaps, mut_count_data, use="complete.obs")
correlations = t(scale(t(correlations)))

write.csv(correlations, "all_corrs.csv")


png(file="/home/mdanb/research/mount_sinai/bing_ren/corrs.png")
h = Heatmap(correlations,
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 6))

dev.off()

