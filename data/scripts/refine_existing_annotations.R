library(optparse)

option_list <- list( 
  make_option("--cores", type="integer")
)

cores = args$cores

###################
# Greenleaf_brain
metadata = read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/GSE162170_atac_cell_metadata_with_cell_names.txt")

glutamatergic_neuron_idx = metadata[, "cell_type"] == "exNeuron"
metadata[glutamatergic_neuron_idx, "refined_cell_type"] = "GluN"

int_neuron_idx = metadata[, "cell_type"] == "intNeuron"
metadata[int_neuron_idx, "refined_cell_type"] = "IN"

unk_idx = metadata[, "cell_type"] == "unknown"
metadata = metadata[!unk_idx, ]

radial_glia_idx = metadata[, "Iterative.LSI.Clusters"] %in% c("c9", "c11")
metadata[radial_glia_idx, "refined_cell_type"] = "RG"

mGPC_idx = metadata[, "Iterative.LSI.Clusters"] == "c10"
metadata[mGPC_idx, "refined_cell_type"] = "mGPC"

IPC_idx = metadata[, "cell_type"] == "IPC"
metadata[IPC_idx, "refined_cell_type"] = "nIPC"

oligo_idx = metadata[, "cell_type"] == "oligo"
metadata[oligo_idx, "refined_cell_type"] = "OPC/Oligo"

mural_idx = metadata[, "cell_type"] == "mural"
metadata[mural_idx, "refined_cell_type"] = "Peric"

microglia_idx = metadata[, "cell_type"] == "microglia"
metadata[microglia_idx, "refined_cell_type"] = "MG"

endothelial_idx = metadata[, "cell_type"] == "endothelial"
metadata[endothelial_idx, "refined_cell_type"] = "EC"

write.table(metadata, 
"/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/Greenleaf_brain_refined_annotations.txt")

###################
# BingRen frontal cortex
metadata = read.csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata/GSE184462_metadata.tsv", 
                      sep="\t",
                      header=T)
metadata = metadata[grepl("Human_brain", metadata[["tissue"]]), ]
astro_1_idx = metadata[, "cell.type"] == "Astrocyte 1"
astro_2_idx = metadata[, "cell.type"] == "Astrocyte 2"
gaba_1_idx = metadata[, "cell.type"] == "GABAergic Neuron 1"
gaba_2_idx = metadata[, "cell.type"] == "GABAergic Neuron 2"
gluta_1_idx = metadata[, "cell.type"] == "Glutaminergic Neuron 1"
gluta_2_idx = metadata[, "cell.type"] == "Glutaminergic Neuron 2"
panc_idx = grepl("Pancreatic", metadata[, "cell.type"])

metadata[astro_1_idx, "cell.type"] = "Astrocyte"
metadata[astro_2_idx, "cell.type"] = "Astrocyte"
metadata[gaba_1_idx, "cell.type"] = "GABAergic Neuron"
metadata[gaba_2_idx, "cell.type"] = "GABAergic Neuron"
metadata[gluta_1_idx, "cell.type"] = "Glutaminergic Neuron"
metadata[gluta_2_idx, "cell.type"] = "Glutaminergic Neuron"
metadata = metadata[!panc_idx, ]

write.table(metadata, 
    "/ahg/regevdata/projects/ICA_Lung/Mohamad/cell_of_origin/data/processed_data/bingren_remove_same_celltype_indexing.tsv")


