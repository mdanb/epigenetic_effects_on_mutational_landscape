library(optparse)
library(readxl)

option_list <- list( 
  make_option("--dataset", type="character"),
  make_option("--annotation", type="character")
)

args = parse_args(OptionParser(option_list=option_list))
# args = parse_args(OptionParser(option_list=option_list), args =
#                  c("--dataset=Greenleaf_colon",
#                    "--annotation=Bingren_remove_same_celltype_indexing"))
# 
dataset = args$dataset
annotation = args$annotation

root = "../metadata"
if (dataset == "Bingren") {
  metadata = read.table("../metadata/GSE184462_metadata.tsv", 
                        sep="\t",
                        header=T)
  if (annotation == "Bingren_remove_same_celltype_indexing") {
    metadata["cell.type"] = gsub(" \\d+", "", metadata[["cell.type"]])
  }
} else if (dataset == "Greenleaf_brain") {
  metadata = read.table("../metadata/GSE162170_atac_cell_metadata.txt.gz", header=1)
  # if (annotation == "Greenleaf_brain_remove_unknown") {
  #   metadata["cell_type"] = gsub(" \\d+", "", metadata[["cell.type"]])
  #   
  #   c("brain (c0$|c1$|c2$|c5|c6|c7|c13|c14)", "brain GluN"),
  #   c("brain (c3|c4|c12|c16|c20)", "brain IN"), 
  #   c("brain c8", "brain nIPC"),
  #   c("brain c9", "brain late RG"),
  #   c("brain c11", "brain early RG"),
  #   c("brain c10", "brain mGPC"),
  #   c("brain c15", "brain OPC/Oligo"),
  #   c("brain c17", "brain Peric"),
  #   c("brain c19", "brain MG"),
  #   c("brain c21", "brain EC"))
  # }
} else if (dataset == "Shendure") {
  metadata = read.table("../metadata/GSE149683_File_S2.Metadata_of_high_quality_cells.txt",
                        sep="\t",
                        header=TRUE)
  if (annotation == "Shendure_remove_unknown_unsure") {
    metadata = metadata[!grepl("Unknown", metadata[["cell_type"]]), ]
    metadata = metadata[!grepl("\\?", metadata[["cell_type"]]), ]
  }
  else if (dataset == "Greenleaf_pbmc_bm") {
    metadata = read.csv("../metadata/greenleaf_pbmc_bm.txt", row.names=1)
    if (annotation == "Greenleaf_pbmc_bm_CD14-mono_CDlike-T_preB+B-B_late+early-no+distinction_Unk-rm") {
      metadata["cell_type"] = gsub("CLP.*", "CLP", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub("CD14.Mono.*", "Monocytes", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub("CD.*", "T", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub("Pre\\.B", 
                                   "B", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub(".*Eryth", 
                                   "Eryth", metadata[["cell_type"]])
      metadata = metadata[!grepl("Unk", metadata[["cell_type"]]), ]
    } else if (annotation == "intermediate_blood_bm_annotation") {
      metadata["cell_type"] = gsub("CLP.*", "CLP", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub("CD14.Mono.*", "Monocytes", 
                                   metadata[["cell_type"]])
      metadata["cell_type"] = gsub("CD.*", "T", 
                                   metadata[["cell_type"]])
      metadata = metadata[!grepl("Unk", metadata[["cell_type"]]), ]
    }
  } else if (dataset == "Yang_kidney") {
    metadata = read_excel("../metadata/41467_2021_27660_MOESM4_ESM.xlsx", 
                               skip=1)
    if (annotation == "Yang_kidney_remove_cell_number_distinctions") {
      metadata["celltype"] = gsub("PT1|PT2|PT3|PT4|PT5|PT6|PT7", "PT", 
                                   metadata[["celltype"]])
      metadata["celltype"] = gsub("CD_ICA1|CD_ICA2", "ICA", 
                                   metadata[["celltype"]])
      metadata["celltype"] = gsub("DCT1|DCT2", "DCT", 
                                   metadata[["celltype"]])
      metadata["celltype"] = gsub("CD_PC1|CD_PC2", "CD_PC", 
                                   metadata[["celltype"]])
    }
  }
}

write.csv(metadata, paste("../metadata",
                          paste(annotation, "metadata.csv", sep="_"), 
                          sep="/"))
