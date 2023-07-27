metadata_epithelial = read.table("../metadata/epithelial_celltypes_atac.tsv", 
                                     header=1)
metadata_stromal = read.table("../metadata/stromal_celltypes_atac.tsv", 
                                  header=1)
metadata_immune = read.table("../metadata/immune_celltypes_atac.tsv", 
                                 header=1)
metadata = rbind(metadata_epithelial, metadata_stromal, metadata_immune)

old_new_names = list(c("A001-C-104-D-R1", "A001-C-104-D_20200214"),
                      c("A001-C-104-D-R2", "A001-C-104-D_20200811"),
                      c("A001-C-124-D-R1", "A001-C-124-D_20200214"),
                      c("A001-C-124-D-R2", "A001-C-124-D_20200702"),
                      c("A001-C-023-D-R1", "A001-C-023-D_20200214"),
                      c("A001-C-023-D-R2", "A001-C-023-D_20200715"),
                      c("A002-C-010-D-R1", "A002-C-010-D_20200310"),
                      c("A015-C-010-D-R2", "A002-C-010-D_20200702"),
                      c("B004-A-004-D-R2", "B004-A-004-D_20200817"))

rename_cells_in_metadata <- function(metadata, old, new) {
  idx = grepl(old, metadata[["Cell"]])
  metadata[idx, "Cell"] = gsub(old, new, 
                               metadata[idx, "Cell"])
  
  return(metadata)
}

for (old_new in old_new_names) {
  metadata[metadata["Sample"] == old_new[[1]], "Sample"] = old_new[[2]]
  metadata = rename_cells_in_metadata(metadata, old_new[[1]], old_new[[2]])
}

write.csv(metadata, "../metadata/greenleaf_colon_metadata.csv")
