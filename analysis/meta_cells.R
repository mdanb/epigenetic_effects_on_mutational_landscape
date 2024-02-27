# Make data.frame of metacolumns to select groups more easily if needed
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/knnGen.R')

metaGroupName='CellType'
metaGroup = as.character(proj@cellColData[,metaGroupName])
barcodes = as.character(proj$cellNames)
metaGroup_df = data.frame(barcode = barcodes, metaGroup = metaGroup)

# KNN = knnGen2 (
#   ArchRProj = archp,
#   reducedDims = 'IterativeLSI',
#   overlapCutoff = 0,
#   cellsToUse = metaGroup_df$barcode,
#   knnIteration = 300,
#   k = 30,
#   min_knn = 2,
#   seed=10
#   )#))

# ###################
# KNN = lapply (unique(as.character(archp@cellColData[,metaGroupName])), function(x) as.list(knnGen2 (
#   ArchRProj = archp,
#   reducedDims = 'IterativeLSI',
#   overlapCutoff = 0,
#   cellsToUse = metaGroup_df$barcode[metaGroup_df$metaGroup == x],
#   knnIteration = 300,
#   k = 30,
#   min_knn = 2,
#   seed=10
#   )))#))





####Filtering step of Archr to Archr2 Filter cell types less than 30 cells)

celltype_keep = names (table (archp@cellColData$predictedGroup_Un)[table (archp@cellColData$predictedGroup_Un) > 30])
archp2 = archp[archp@cellColData$predictedGroup_Un %in% celltype_keep]
####Again KNN


KNN = lapply (unique(as.character(proj@cellColData[,metaGroupName])), function(x) as.list(knnGen2 (
  ArchRProj = archp2,
  reducedDims = 'IterativeLSI',
  overlapCutoff = 0,
  cellsToUse = metaGroup_df$barcode[metaGroup_df$metaGroup == x],
  knnIteration = 300,
  k = 30,
  min_knn = 2,
  seed=10
)))#))

KNN = unlist(KNN, recursive = F)