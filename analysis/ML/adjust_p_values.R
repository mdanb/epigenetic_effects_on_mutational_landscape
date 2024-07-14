
first_second = read.csv("models/XGB/p_values_feat_imp_5_feats.csv", row.names=1)
colnames(first_second) = gsub("\\.","-", colnames(first_second))

first_third = read.csv("models/XGB/p_values_feat_imp_third_5_feats.csv", row.names=1)
colnames(first_third) = gsub("\\.","-", colnames(first_third))

lusc=c(first_second["Basal, Lung D5", "Lung-SCC"], 
       first_third["Basal, Lung D5", "Lung-SCC"])
p.adjust(lusc, method="BH")

meso=c(first_second["Mesothelium,\nLung D5", "epithelioid_waddell"], 
       first_third["Mesothelium,\nLung D5", "epithelioid_waddell"])
p.adjust(meso, method="BH")

sclc=c(first_second["Basal, Lung D5", "SCLC"], 
       first_third["Basal, Lung D5", "SCLC"])
p.adjust(sclc, method="BH")


first_second_top_feat = read.csv("models/XGB/p_values_top_feature.csv", row.names=1)
colnames(first_second_top_feat) = gsub("\\.","-", colnames(first_second_top_feat))
first_third_top_feat = read.csv("models/XGB/p_values_top_feature_third_top.csv", row.names=1)
colnames(first_third_top_feat) = gsub("\\.","-", colnames(first_third_top_feat))

aml = c(first_second_top_feat["bonemarrow GMP GL_BlBm", "Myeloid-AML"], 
       first_third_top_feat["bonemarrow GMP GL_BlBm", "Myeloid-AML"])
p.adjust(aml, method="BH")
