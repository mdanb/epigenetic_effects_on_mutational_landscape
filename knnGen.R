#Run This to accesss hidden functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}


knnGen2 = function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  corCutOff = 0.75,
  dimsToUse = 1:30,
  knnIteration = 500,
  k = 100, # automatically set if group is specified
  min_knn = 2,
  overlapCutoff = 0,
  seed = 1,
  cellsToUse = NULL
  )
{
  set.seed (seed)
  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix

  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,, drop=F]

  # if number of k doesnt meet minimum
  if (sum(keepKnn == 0) < min_knn)
    {
    slice = function(x,n) {
      N = length(x);
      lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
     }
    message ('KNN search failed. Group cells by k')  
    knnObj =  sapply (slice (rownames(rD), k), function(x) match (x,rownames(rD)))  
    knnObj = knnObj[sapply(knnObj, length) == k]
    knnObj = do.call (rbind, knnObj)
    }

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
 return (knnObj)
}

### KNN generator function ###
# k # number of neighbors per knnIteration NOTE: this get overwritten if group is specified
# knnIteration # number of cells sampled from which k neighbors are searched NOTE: ignored if gourp is specified
# min.cells_in_group # ignored if group is NULL
# min_knn_cluster # ignored if group is NULL
knnGen = function (
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  corCutOff = 0.75,
  dimsToUse = 1:30,
  knnIteration = 500,
  k = 100, # automatically set if group is specified
  overlapCutoff = 0.8,
  seed = 1,
  cellsToUse = NULL,
  group = NULL,
  min.cells_in_group = 50,
  min_knn_cluster = 2,
  max_attempts = 1000
  )
  {

  require (progress)
  set.seed (seed)

  slice = function(x,n) {
    N = length(x);
    lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
  }
  #Get Reduced Dims
  rD = getReducedDims (ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  
  if(!is.null(cellsToUse))
  {
  rD = rD[cellsToUse, ,drop=FALSE]
  rD = as.data.frame (rD)
  ArchRProj = ArchRProj[cellsToUse,]  
  }
  
  if (!is.null (group))
    {
    meta_group = getCellColData (ArchRProj)[,group]  
    meta_group = table (meta_group[!meta_group %in% names(table (meta_group)[table(meta_group) <= min.cells_in_group])])
    #k = floor (min (meta_group) / min_knn_cluster)
    #message (paste('k set to', k, 'to allow at least',min_knn_cluster,'knn group in each cluster'))
    knnObj_l = list()
    for (i in names(meta_group))
      {
      keepKnn = NULL  
      l = 1
      message (paste('find at least',min_knn_cluster,'knn in cluster',i))
      while (sum(keepKnn == 0) < min_knn_cluster)
        {
        #Subsample
        idx <- unique (sample (seq_len(meta_group[i]), knnIteration, replace = !meta_group[i] >= knnIteration))
        #KNN Matrix
        rDs = rD[rownames(rD) %in% ArchRProj$cellNames[as.character(ArchRProj@cellColData[,group]) == i],]
        knnObj <- .computeKNN (data = rDs, query = rDs[idx,], k = k)
        #Determin Overlap
        keepKnn <- determineOverlapCpp (knnObj, floor(overlapCutoff * k))

        #Keep Above Cutoff
        knnObj_group = knnObj[keepKnn==0,, drop=F] # remove knn with too many overlapping cells
        l = l + 1
        if (max_attempts == l) 
          {
          message ('KNN search failed. Partition the cell group')  
          knnObj_group = lapply (slice (rownames(rDs), k), function(x) match (x,rownames(rDs)))  
          knnObj_group = knnObj_group[sapply(knnObj_group, length) == k]
          knnObj_group = do.call (rbind, knnObj_group)
          break
          }
        }
      knnObj_l[[i]] = apply (knnObj_group, 2, function(x) rownames(rDs)[x])
      if (!is.matrix(knnObj_l[[i]])) knnObj_l[[i]] = t(as.matrix (knnObj_l[[i]]) )
      rownames (knnObj_l[[i]]) = paste0(i, 'KNN',1:nrow(knnObj_l[[i]]))
      }    
  knnObj = Reduce (rbind, knnObj_l)
  knn_names = rownames (knnObj) 
  knnObj = lapply (seq_len(nrow(knnObj)), function(x) knnObj[x,]) %>% SimpleList
  names (knnObj) = knn_names
  
  } else {
  message ('Run cluster-agnostic knn')  
  #Subsample
  idx <- sample (seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  knnObj <- .computeKNN (data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  keepKnn <- determineOverlapCpp (knnObj, floor(overlapCutoff * k))

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,] # remove knn with too many overlapping cells
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
  rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  names (knnObj) = paste0 ('KNN_',seq_along (knnObj))
  }
  return (knnObj)
}








