

metadata_epithelial = read.table("../metadata/epithelial_celltypes_atac.tsv", 
                                     header=1)
metadata_stromal = read.table("../metadata/stromal_celltypes_atac.tsv", 
                                  header=1)
metadata_immune = read.table("../metadata/immune_celltypes_atac.tsv", 
                                 header=1)
metadata = rbind(metadata_epithelial, metadata_stromal, metadata_immune)
metadata["Sample"] = gsub("-R[1-2]$", "", metadata[["Sample"]])

write.csv(metadata, "../metadata/greenleaf_colon_metadata.csv")
