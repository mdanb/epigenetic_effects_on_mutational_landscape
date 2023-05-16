import pandas as pd

df_1 = pd.read_csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/"
                   "mutation_data/mesothelioma_WGS_Waddell.csv", index_col=0)

df_2 = pd.read_csv("/broad/hptmp/bgiotti/BingRen_scATAC_atlas/raw_dir/"
                   "mutation_data/mesothelioma_p786_p848.csv", index_col=0)

df = df_1.join(df_2)
df.to_csv("processed_data/meso_mut_count_data.csv")