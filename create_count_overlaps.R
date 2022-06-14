library(tools)                                                                  
library(rtracklayer)                                                            
library(stringr)                                                                
library(data.table)                                                             
library(dplyr)                                                                  
library(tibble)                                                                 

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied", call.=FALSE)
} else {
  CELL_NUMBER_FILTER = args[1]
} 

load('/ahg/regevdata/projects/ICA_Lung/Bruno/celloforigin_prj/hg19.1Mb.ranges.Polak.Nature2015.RData')
#load('hg19.1Mb.ranges.Polak.Nature2015.RData')

interval.ranges                                 

dir.create("count_overlap_data")                                                

metadata = read.table("raw_dir/metadata/GSE184462_metadata.tsv", sep="\t", header=TRUE)  

for (file in setdiff(list.files("raw_dir/bed_files/"), list.dirs("raw_dir/bed_files", 
		     recursive = FALSE, full.names = FALSE))) {                                    
  filename = paste("count_filter", CELL_NUMBER_FILTER,  
                   "count_overlaps", paste(file_path_sans_ext(file, TRUE),
	                 "rds", sep="."), sep="_")
	if (!file.exists(filename)) {
		print(paste("Processing", file, sep= " "))
		sample = import(paste("raw_dir", "bed_files", file, sep="/"), format="bed")
		if (grepl("frontal_cortex", file)) {
			sample_name = "frontal_cortex"
		} else {
			sample_name = str_extract(file, pattern="_.*_SM")              
			sample_name = substr(sample_name, 2, nchar(sample_name) - 3)    
                }
                filtered_metadata = metadata[metadata$sample %like% sample_name, ]
                sample_cellIDs = filtered_metadata[["cellID"]]                  
                sample_barcodes_in_metadata = unlist(lapply(sample_cellIDs, function(x) unlist(strsplit(x, '\\+'))[2]))
                sample = sample[sample$name %in% sample_barcodes_in_metadata]   
                sample_idx_in_metadata = match(sample$name, sample_barcodes_in_metadata)
                sample$cell_type = filtered_metadata[unlist(sample_idx_in_metadata), "cell.type"]
                sample = as.tibble(sample)                                 
                counts_per_cell_type = sample %>%                               
                                       group_by(cell_type) %>%                 
                                       summarize(n_cells = n_distinct(name))    

                cell_type_keep = counts_per_cell_type[counts_per_cell_type$n_cells >= CELL_NUMBER_FILTER, ]$cell_type
                sample = sample %>%                                             
                         filter(cell_type %in% cell_type_keep)                  

                grl_in = sample %>%                                             
		         group_split(cell_type) %>%                             
			 lapply(makeGRangesFromDataFrame, keep.extra.columns=TRUE)
		grl = GRangesList(grl_in)
		count_overlaps = lapply(grl, function(x) countOverlaps (interval.ranges, x))
		cell_types = unlist(lapply(grl, function(x) unique(x$cell_type)))
		names(count_overlaps) = cell_types
		saveRDS(count_overlaps, paste("count_overlap_data", filename, sep="/"))
	}
}
