import scanpy

df = scanpy.read_h5ad("../211201Fetal_lung_ATAC.h5ad")
cell_types = df.obs[["cell_type"]]
cell_types.to_csv("rawlins_fetal_lung_metadata.csv")
