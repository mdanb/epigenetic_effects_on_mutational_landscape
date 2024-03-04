library(Seurat)
library(ggplot2)   
library(ggpubr)    

### Cell cycle analysis ###
# Generate barplots or boxplots of meta groups proportions specified
# meta_groups: vector or meta group names:
# 1) first meta_group is the group on which proportions are calculated
# 2) second meta_group split first meta_group on x axes  
#   3) third meta_group will group barplots separately
# if splits include only one value runs barplot instead
cellComp = function (
    seurat_obj = NULL, 
    metaGroups = NULL, # vector of at least 3 metaGroups e.g. c('orig.ident','celltypes','celltypes'),
    plot_as = 'box', # box or bar 
    pal = NULL,
    prop = TRUE,
    ptable_factor = 1, # specify which column of the data.frame or seurat object metadata should be used to compute proportions
    facet_ncol = 20,
    facet_scales = 'free',
    subset_prop = NULL, # subset prop table by any group in any column
    removeNA = TRUE,
    returnDF = FALSE
)
{
  if (is.data.frame (seurat_obj))
  {
    meta_groups_df = seurat_obj[,metaGroups]    
  } else {
    meta_groups_df = seurat_obj@meta.data[,metaGroups]
  }
  # Refactor to remove 0 groups
  #meta_groups_df =  as.data.frame(lapply(unclass(meta_groups_df),as.character),stringsAsFactors=T)
  if(is.null(pal)) pal = rainbow (length(unique(meta_groups_df[,2])))
  #if(is.null(pal) & plot_as == 'box') pal = rainbow (length(unique(meta_groups_df[,3])))
  if (prop)
  {
    ccomp_df = as.data.frame (prop.table (table (meta_groups_df),ptable_factor))
    ccomp_df = na.omit (ccomp_df) # this is to remove NaN somehow produced from the line above 
  } else {
    ccomp_df = as.data.frame (table (meta_groups_df))   
  }
  
  if(removeNA) ccomp_df = ccomp_df[ccomp_df$Freq != 0, ] # remove 0s from proportions
  if (!is.null (subset_prop)) 
  {
    subset_col = unlist(sapply (seq(ncol(ccomp_df)), function(x) if(any(ccomp_df[,x] %in% subset_prop)) colnames(ccomp_df)[x]))
    ccomp_df = ccomp_df[ccomp_df[,subset_col] %in% subset_prop,]
  }
  #colnames (ccomp_df) = c(paste0('Var_',seq_along(metaGroups)), 'proportion')  
  if (plot_as == 'box')
  {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[2], y= 'Freq')) +
      theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups > 2)) p = p + geom_boxplot(aes_string (fill= metaGroups[3]), outlier.size=.2, alpha = 0.7, lwd=.2) 
    else p = p + geom_boxplot(aes_string (fill= metaGroups[2]), outlier.size=.2, alpha = 0.7, lwd=.2)       
    if (length(metaGroups) > 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[4])), scales=facet_scales, ncol=facet_ncol)
  }   
  if (plot_as == 'bar')
  {
    p = ggplot (ccomp_df, aes_string (x= metaGroups[1], y= 'Freq')) +
      geom_bar(position="stack", stat="identity", aes_string(fill= metaGroups[2])) +
      #geom_bar(position="dodge", stat="identity", aes_string(fill= metaGroups[2])) +
      theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_fill_manual (values= pal) + xlab (metaGroups[2]) + ylab (ifelse (prop, 'proportion','counts'))
    if (length(metaGroups) == 3) p = p + facet_wrap (as.formula(paste("~", metaGroups[3])), scales=facet_scales, ncol=facet_ncol)
  }
  if (returnDF) return(ccomp_df) else 
    return (p)
}

### Variables needed ###
metaGroupName1 = "sample_ID"
metaGroupName2 = "Cluster"

# cc_transfer # assign cell identity of cycling cells using label transfer. Logical
srt = readRDS("../data/scRNAseq/tosti_pancreas/adult_pancreas_2020.rds")
if (!'cc' %in% colnames (srt@meta.data)) {
  cc.genes <- readLines('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/regev_lab_cell_cycle_genes.txt')
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]
  
  message ("Calculate cell cycle score with Seurat CellCycleScoring function")
  srt = CellCycleScoring(object = srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srt$cc = srt$S.Score + srt$G2M.Score
} else {
  message ("Cell cycle scores found...generating plots")   
}  

srt$Phase2 = ifelse (srt$Phase == 'G1','non-cycling','cycling')
metaGroupNames = c(metaGroupName2, 'Phase2')
# ccc_box1 = cellComp (
#   seurat_obj = srt, 
#   metaGroups = metaGroupNames,
#   plot_as = 'box',
# #  pal = pal1,
#   prop = TRUE,
#   ptable_factor = c(1,3),
#   subset_prop = 'cycling',
#   facet_ncol = 6,
#   facet_scales = 'free'
#   ) + theme_classic()#+ NoLegend() + ggtitle (i)

ccc_bar1 = cellComp(
  seurat_obj = srt, 
  metaGroups = metaGroupNames,
  plot_as = 'bar',
  pal = c(cycling = 'blue', `non-cycling`= 'grey'),
  prop = TRUE,
  ptable_factor = c(1),
  #subset_prop = 'cycling',
  facet_ncol = 6,
  facet_scales = 'free'
) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=25),
        axis.text.y = element_text(size=25))
#+ NoLegend() + ggtitle (i)
# fn = paste(metaGroupNames, collapse='_'), '_composition.pdf')
pdf(paste0("../figures/", paste(metaGroupNames, collapse='_'),
            '_composition.pdf'), width=20, height=15)
print ((ccc_bar1))# + plot_layout(widths=c(1,1))
dev.off()

p = FeaturePlot(srt, features="G2M.Score")
pdf("../figures/temp.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()

p = DimPlot(srt, group.by="Phase")
pdf("../figures/phase.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()


p = DimPlot(srt, group.by="Cluster")
pdf("../figures/cluster.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()

p = FeaturePlot(srt, features="cc")
pdf("../figures/cc.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()

"TOP2A"
"MKI67"

p = FeaturePlot(srt, features="TOP2A")
pdf("../figures/top2a.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()

p = FeaturePlot(srt, features="MKI67")
pdf("../figures/mki67.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()

p = FeaturePlot (srt, 'TOP2A', combine=FALSE) + 
  scale_colour_gradientn (colours=viridis::turbo(10))
pdf("../figures/top2a.pdf", width=20, height=15)
print ((p))# + plot_layout(widths=c(1,1))
dev.off()
