library(tools)
library(rtracklayer)

load('raw_dir/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData')

dir.create("polak_count_overlap_data")

for (file in list.files("polak_raw_data/")) {
	filename = paste("count_overlaps", paste(file_path_sans_ext(file, TRUE), 
			 "rds", sep="."), sep="_")
	if (!file.exists(filename)) {
		tryCatch({	
			bed_file = import(paste("polak_raw_data", file, sep="/"), format="bed")
			count_overlaps = countOverlaps(interval.ranges, bed_file)
			saveRDS(count_overlaps, paste("polak_count_overlap_data", filename, sep="/"))
			unlink(paste("polak_raw_data", file, sep="/"))
		},
		error = function(err) {
		})
	}
	else {
		unlink(paste("polak_raw_data", file, sep="/"))
	}
}
