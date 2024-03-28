library(optparse)
library(gridExtra)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)

source("ML_utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("--datasets"), type="character")
parser <- add_option(parser, c("--cancer_types"), type="character")
parser <- add_option(parser, c("--cell_number_filter"), type="integer")
parser <- add_option(parser, c("--tss_fragment_filter"),
                     type="character", default="-1")
parser <- add_option(parser, c("--tissues_to_consider"), 
                     type="character", default="all")
parser <- add_option(parser, c("--ML_model"), type="character")
parser <- add_option(parser, c("--annotation"), 
                     type="character", default="default_annotation")
parser <- add_option(parser, c("--feature_importance_method"), type="character")
parser <- add_option(parser, c("--folds_for_test_set"), type="character")
parser <- add_option(parser, c("--plot_bin_counts"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--plot_bins_volcano"), action="store_true", 
                     default=F)
parser <- add_option(parser, c("--hundred_kb"), action="store_true", 
                     default=F)

args = parse_args(parser)
# 
args = parse_args(parser, args =
                    c("--datasets=Bingren",
                      "--cancer_types=Skin-Melanoma",
                      "--cell_number_filter=100",
                      "--ML_model=XGB",
                      "--annotation=Bingren_remove_same_celltype_indexing",
                      "--feature_importance_method=permutation_importance",
                      "--folds_for_test_set=1-10",
                      "--hundred_kb",
                      "--plot_bins_volcano"))

args = parse_args(parser, args =
                    c("--datasets=Tsankov",
                      "--cancer_types=Lung-SCC",
                      "--cell_number_filter=100",
                      "--ML_model=XGB",
                      "--annotation=test_annotation",
                      "--feature_importance_method=permutation_importance",
                      "--folds_for_test_set=1-10"))
cancer_types = args$cancer_types
cancer_types = unlist(strsplit(cancer_types, split = ","))
datasets = unlist(strsplit(args$datasets, split = ","))
datasets = sort(datasets)
cell_number_filter = args$cell_number_filter
tss_fragment_filter = unlist(strsplit(args$tss_fragment_filter, split = ","))
annotation = args$annotation
tissues_to_consider = strsplit(args$tissues_to_consider,  split=",")
ML_model = args$ML_model
feature_importance_method = args$feature_importance_method

folds_for_test_set = args$folds_for_test_set
folds_for_test_set = unlist(strsplit(args$folds_for_test_set, split = "-"))
folds_for_test_set = seq(folds_for_test_set[1], folds_for_test_set[2])
plot_bin_counts = args$plot_bin_counts
hundred_kb = args$hundred_kb
plot_bins_volcano = args$plot_bins_volcano

plot_dist <- function(X, bin_names) {
  df <- data.frame(index = 1:nrow(X), 
                          value = X)
  use_names = seq(1, length(df$index), by=10)
  use_names = append(use_names, nrow(X))
  p = ggplot(df) +
    geom_bar(aes(x=index, y=value),
             stat="identity") +
    xlab("Bin") +
    ylab("Counts") +
    geom_text(x=nrow(df) / 2, y=max(df[["value"]]), 
              aes(label=paste0("n=", sum(df[["value"]]))),
              size=10) +
    scale_x_continuous(name="Row Names", breaks=df$index[use_names], 
                       labels=bin_names[use_names]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}


# scatac_counts_plots <- list()
# mut_counts_plots <- list()

for (cancer_type in cancer_types) {
  scATAC_sources = construct_sources_string(datasets)
  if (plot_bins_volcano) {
    dir = construct_dir(scATAC_sources,
                        cell_number_filter,
                        tss_fragment_filter,
                        annotation,
                        "all_seeds",
                        -1,
                        ML_model,
                        cancer_type,
                        hundred_kb)
    if (hundred_kb) {
      load("../../data/100kb_interval_ranges.Rdata")
    }
    else {
      load("../../data/mutation_data/hg19.1Mb.ranges.Polak.Nature2015.RData")
    }
    
    fp = paste("../../figures", dir, "errors_df.csv", sep="/")
    df = read.csv(fp)
    colnames(df)[1] = "names"
    chr_to_range = data.frame(interval.ranges@ranges)
    df = left_join(df, chr_to_range)
    # df["neglog_q"] = -log10(df["absolute_error_q.value"])
    # df = df %>% arrange(desc(percent_error))
      
    # ggplot(df) + 
    #   geom_histogram(aes(x=absolute_error))
    # 
    # ggplot(df) + 
    #   geom_histogram(aes(x=log(abs(percent_error)+1)))
    # 
    # ggplot(df) + 
    #   geom_histogram(aes(x=percent_error))
    
    # df = df %>%
    #       filter(percent_error_rejected == "True")
    
    # df = df %>% arrange(desc(normalized_absolute_error))
    # df = head(df, 1000)
    # ggplot(df) + 
    #   geom_point(aes(x=normalized_absolute_error, y=neglog_q)) +
    #   ylab("Mutation Enrichment, -log10 FDR")
    # 
    # df[, c("names", "start", "end", "neglog_q")]
  }
  
  chr = unlist(lapply(strsplit(df$names, "\\."), "[", 1))
  df["chr"] = chr
  bins = GRanges(seqnames=df$chr, IRanges(start=df$start, end=df$end),
                 avg_absolute_error=df$avg_absolute_error,
                 median_absolute_error=df$median_absolute_error,
                 avg_pred=df$avg_prediction,
                 median_pred=df$median_prediction,
                 actual=df$actual)
  bins_df = as.data.frame(bins)
 
  ggplot(bins_df) +
    geom_histogram(aes(x=actual))
  
  ggplot(bins_df) +
    geom_histogram(aes(x=avg_pred))
  
  ggplot(bins_df) +
    geom_histogram(aes(x=median_pred))
  
  bins_df = bins_df %>% 
              mutate(avg_percent_error = 100*(actual-avg_pred)/(actual+1),
                     dinf_avg = 100*(actual-avg_pred) / max(abs(actual),
                                                      abs(avg_pred)),
                     d1_avg = 200*(actual-avg_pred) / (abs(actual) + abs(avg_pred)),
                     dg_avg = 100*(actual-avg_pred) / (sqrt(abs(actual*avg_pred)+1))
                     )
  print(max(bins_df["avg_percent_error"]))
  print(min(bins_df["avg_percent_error"]))
  print(max(bins_df["dinf_avg"]))
  print(min(bins_df["dinf_avg"]))
  print(max(bins_df["dg_avg"]))
  print(min(bins_df["dg_avg"]))
  print(max(bins_df["d1_avg"]))
  print(min(bins_df["d1_avg"]))
  
  ggplot(bins_df) +
    geom_histogram(aes(x=avg_percent_error)) 

  ggplot(bins_df) +
    geom_histogram(aes(x=dg_avg)) +
    xlim(-750, 300)
  
  ggplot(bins_df) +
    geom_histogram(aes(x=dinf_avg))

  ggplot(bins_df) +
    geom_histogram(aes(x=d1_avg))
  
  ggplot(bins_df) +
    geom_point(aes(x=avg_absolute_error, y=avg_percent_error, 
                   size=actual, color=avg_pred)) +
    xlab("Average Absolute Error") +
    ylab("Average Percent Error") 
    # ylim(-750, 300)
  
  ggplot(bins_df) +
    geom_point(aes(x=avg_absolute_error, y=dg_avg, size=actual, color=avg_pred)) +
    xlab("Average Absolute Error") +
    ylab("Average dg") #+
   # ylim(-750, 300)
  
  ggplot(bins_df) +
    geom_point(aes(x=avg_absolute_error, y=dinf_avg), size=0.2) +
    xlab("Average Absolute Error") +
    ylab("Average dg") 
  
  ggplot(bins_df) +
    geom_point(aes(x=avg_absolute_error, y=dinf_avg, size=actual, color=avg_pred)) +
    xlab("Average Absolute Error") +
    ylab("Average dg") 
  
  ggplot(bins_df) +
    geom_point(aes(x=avg_absolute_error, y=d1_avg,size=actual, color=avg_pred))+
    xlab("Average Absolute Error") +
    ylab("Average d1") 
  # basal_se = import.bed("../../data/40_hubscut4tss2000min3/AT2.bed") 
  # seqlevels(basal_se) = paste0("chr", levels(seqnames(basal_se)))
  # basal_se_total = findOverlaps(basal_se, significant_ranges)
  # total_matches = length(basal_se_total)
  # total_num_bins = 21280
  # df = head(df, 1000)
  # significant_ranges = GRanges(seqnames=df$chr, IRanges(start=df$start, end=df$end))
  # basal_se_top_1000 = findOverlaps(basal_se, significant_ranges)
  # top_1000_num_matches = length(basal_se_top_1000)
  
  # chr_ranges = read.csv("../../data/processed_data/chr_ranges.csv")
  chr_keep_100kb = read.csv("../../data/processed_data/chr_keep_100kb.csv",
                            row.names = 1)
  load('../../data/100kb_interval_ranges.Rdata')
  
  cosmic_genes = read.csv("../../data/cosmic_genes.csv")
  gr = strsplit(cosmic_genes$Genome.Location, split=":")
  cosmic_genes = cosmic_genes[lapply(gr, "[", 2) != "-", ]
  gr = strsplit(cosmic_genes$Genome.Location, split=":")
  chr = lapply(gr, "[", 1)
  chr = paste0("chr", chr)
  start_end = strsplit(unlist(lapply(gr, "[", 2)), split="-")
  start = as.numeric(unlist(lapply(start_end, "[", 1)))
  end = as.numeric(unlist(lapply(start_end, "[", 2)))
  cosmic_gene_ranges = GRanges(seqnames=chr, IRanges(start=start, end=end))
  overlaps = findOverlaps(significant_ranges, cosmic_gene_ranges)
  overlaps = as.data.frame(overlaps)
  # hits_ours = df[overlaps$queryHits, ]
  # significant_ranges
  hits_cosmic = cosmic_genes[overlaps$subjectHits, ]
  # hits_cosmic["percent_error"] = significant_ranges[overlaps$queryHits]$percent_error
  # hits_cosmic["absolute_error"] = significant_ranges[overlaps$queryHits]$absolute_error
  # hits_cosmic = hits_cosmic %>% 
  #                 group_by(Gene.Symbol) %>% 
  #                 summarise(avg_percent_error = mean(percent_error),
  #                           avg_absolute_error = mean(absolute_error))
  # hits_cosmic = hits_cosmic %>%
  #                 left_join(cosmic_genes) %>% 
  #                 arrange(desc(avg_percent_error))
  # 
  # N = nrow(hits_cosmic)
  # cancer_name_cosmic = "lung$|lung,|lung cancer|lung SCC"
  # cancer_name_cosmic = "lung SCC"
  
  # K = sum(unlist(lapply(hits_cosmic["Tumour.Types.Somatic."], 
  #            function(x) { grepl(cancer_name_cosmic, x, ignore.case = T) })))
  # n = floor(0.1 * N)
  # k = sum(unlist(lapply(hits_cosmic[1:n, "Tumour.Types.Somatic."], 
  #            function(x) { grepl(cancer_name_cosmic, x, ignore.case = T) })))
  # p_value <- phyper(k, K, N-K, n, lower.tail = F)
                            
  # a = distinct(hits_cosmic)
  # write.csv(x=df, file=paste(cancer_type, "all_hits_ours.csv", sep="_"))
  # write.csv(x=hits_ours, file=paste(cancer_type, "hits_ours.csv", sep="_"))
  # write.csv(x=hits_cosmic, file=paste(cancer_type, "hits_cosmic_ours.csv", sep="_"))
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  broads <- GenomicFeatures::genes(txdb, columns="gene_id")
  
  ###
  # tp53 = dplyr::filter(broads, gene_id == "7157")
  # chr17 = significant_ranges %>% filter(seqnames == "chr17")
  # 
  # as.data.frame(chr17)[start(ranges(chr17)) > 7500000 & 
  #                        end(ranges(chr17)) < 7600000 ,]
  # findOverlaps(significant_ranges, tp53)

  
  ###
  
  overlaps = as.data.frame(findOverlaps(bins, broads))
  entrez_ids <- broads[overlaps$subjectHits]$gene_id
  errors = as.data.frame(bins[overlaps$queryHits])[c("avg_percent_error",
                                                     "avg_absolute_error",
                                                     "avg_pred",
                                                     "actual")]
  
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys=entrez_ids,
                         column="SYMBOL",
                         keytype="ENTREZID",
                         multiVals="first")
  broads = as.data.frame(broads)
  # rownames(broads) = broads$gene_id
  
  # if (cancer_type == "Lung-AdenoCA") {
  #   gene_symbols[is.na(gene_symbols)] = "FBXO22-AS1 (deprecated)"
  # }
  
  genes = data.frame(symbol=gene_symbols)
  genes = cbind(genes, broads[entrez_ids, ], errors)
  genes = genes %>%
            # filter(!is.na(symbol)) %>%
            group_by(gene_id, symbol) %>% 
            summarise(avg_percent_error = mean(avg_percent_error),
                      avg_absolute_error = mean(avg_absolute_error))
  
  # genes = genes[!duplicated(cbind(genes$symbol, genes$entrez_id)), ]
  # genes = cbind(genes, broads[genes[["entrez_id"]], ])
  melanoma_genes = unique(hits_cosmic[grepl("melanoma",
                    hits_cosmic[["Tumour.Types.Somatic."]]),][["Gene.Symbol"]])
  melanoma_entrez_ids = mapIds(org.Hs.eg.db,
                                   keys=melanoma_genes,
                                   column="ENTREZID",
                                   keytype="SYMBOL",
                                   multiVals="first")
  
  
  lung_genes = read.csv("../../data/hallmark_genes.csv")
  
  lung_adenoca_entrez_ids = mapIds(org.Hs.eg.db,
                               keys=lung_genes[["Lung.AdenoCA"]],
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
  
  # FTSJD1: CMTR2
  # KARS: KARS1
  # MLL3: KMT2C
  
  lung_genes[lung_genes[["Lung.AdenoCA"]] == "KARS", "Lung.AdenoCA"] = "KARS1"
  lung_genes[lung_genes[["Lung.AdenoCA"]] == "FTSJD1", "Lung.AdenoCA"] = "CMTR2"
  lung_genes[lung_genes[["Lung.AdenoCA"]] == "MLL3", "Lung.AdenoCA"] = "KMT2C"
  lung_adenoca_entrez_ids = mapIds(org.Hs.eg.db,
                                   keys=lung_genes[["Lung.AdenoCA"]],
                                   column="ENTREZID",
                                   keytype="SYMBOL",
                                   multiVals="first")
  lung_adenoca_entrez_ids = lung_adenoca_entrez_ids[lung_adenoca_entrez_ids %in% 
                                                      genes[["gene_id"]]]
  
  lung_scc_entrez_ids = mapIds(org.Hs.eg.db,
                               keys=lung_genes[["Lung.SCC"]][lung_genes[["Lung.SCC"]] != ""],
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")

  # MLL2: KMT2D
  lung_genes[lung_genes[["Lung.SCC"]] == "MLL2", "Lung.SCC"] = "KMT2D"
  lung_scc_entrez_ids = mapIds(org.Hs.eg.db,
                               keys=lung_genes[["Lung.SCC"]][lung_genes[["Lung.SCC"]] != ""],
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")
  
  lung_scc_entrez_ids = lung_scc_entrez_ids[lung_scc_entrez_ids %in% 
                                              genes[["gene_id"]]]
  
  conduct_hypergeom_test <- function(underestimated, rank_by, all_genes, 
                                     gene_subset, top_a_perc) {
    genes = all_genes %>%
              filter(!is.na(symbol))
    if (underestimated) {
        genes = genes %>% 
                  arrange(desc(!!sym(rank_by)))
    } else {
        genes = genes %>% 
            arrange(!!sym(rank_by))
    }
    
    N = nrow(genes)
    K = length(gene_subset)
    n = floor(top_a_perc * N)
    k = sum(genes[["gene_id"]][1:n] %in% gene_subset)
    p_value <- phyper(k, K, N-K, n, lower.tail = F)
    print(paste(paste0(k,"/",K), "genes in top", paste0(top_a_perc * 100, "%")))
    print(paste0("p=", p_value))
    top_n = genes[1:n, ][, c("gene_id", "symbol")]
    print("genes:")
    print( top_n[["symbol"]][top_n[["gene_id"]] %in% gene_subset])
  }
  
  gene_subset = melanoma_entrez_ids
  gene_subset = lung_scc_entrez_ids
  gene_subset = lung_adenoca_entrez_ids
  # rank_by = "avg_absolute_error"
  rank_by = "avg_percent_error"

  conduct_hypergeom_test(underestimated=F, rank_by=rank_by, 
                         all_genes=genes, gene_subset=gene_subset,
                         top_a_perc = 0.1)
  
  significant_ranges = sort(significant_ranges, by=~percent_error)
  percent = 5
  n = floor(5/100*length(significant_ranges))
  top_n_bins = significant_ranges[1:n]
  temp = length(significant_ranges)-n+1
  bottom_n_bins = significant_ranges[temp:length(significant_ranges)]
  strand = "minus"
  fn_top = paste(cancer_type, "top" , percent, "percent_bins", "strand",
                 paste0(strand, ".bed"), sep="_")
  fp_top = paste("post_modeling_results", fn_top, sep="/")
  fn_bottom = paste(cancer_type, "bottom" , percent, "percent_bins", 
                    "strand", paste0(strand, ".bed"), sep="_")
  fp_bottom = paste("post_modeling_results", fn_bottom, sep="/")
  
  if (strand == "plus") {
    strand = "+"
  } else if (strand == "minus") {
    strand = "-"
  }
  strand(top_n_bins) = strand
  strand(bottom_n_bins) = strand
  export(top_n_bins, fp_top, format="bed")
  export(bottom_n_bins, fp_bottom, format="bed")
  # N = nrow(genes)
  # K = length(lung_scc_entrez_ids)
  # n = floor(0.1 * N)
  # k = sum(genes[["gene_id"]][1:n] %in% lung_scc_entrez_ids)
  # p_value <- phyper(k, K, N-K, n, lower.tail = F)
  # 
  # K = sum(unlist(lapply(hits_cosmic["Tumour.Types.Somatic."], 
  #                       function(x) { grepl(cancer_name_cosmic, x, ignore.case = T) })))
  # n = floor(0.1 * N)
  # k = sum(unlist(lapply(hits_cosmic[1:n, "Tumour.Types.Somatic."], 
  #                       function(x) { grepl(cancer_name_cosmic, x, ignore.case = T) })))
  # p_value <- phyper(k, K, N-K, n, lower.tail = F)
  
  # write.csv(x=genes, file=paste(cancer_type, "all_gene_hits.csv", sep="_"),
  #           row.names = F)
  
  
  if (plot_bin_counts) {
    for (fold in folds_for_test_set) {
    # for (tss_filter in tss_fragment_filter) {
      dir = construct_dir(scATAC_sources,
                    cell_number_filter,
                    tss_fragment_filter,
                    annotation,
                    1,
                    fold,
                    ML_model,
                    cancer_type, 
                    hundred_kb)

      X_test = read.csv(paste(dir, "X_test.csv", sep="/"), row.names = 1,
                        check.names = FALSE)
      y_test = read.csv(paste(dir, "y_test.csv", sep="/"), row.names = 1)
      backwards_elim_dir = construct_backwards_elim_dir(cancer_type, 
                                                        scATAC_sources, 
                                                        cell_number_filter,
                                                        tss_fragment_filter, 
                                                        annotation,
                                                        tissues_to_consider, 
                                                        ML_model,
                                                        1,
                                                        hundred_kb = hundred_kb,
                                                        fold_for_test_set=fold,
                                                        test=T)
      fp = paste(backwards_elim_dir, 
              "top_features_iteration_19_by_permutation_importance.txt", sep="/")
      con = file(fp, "r")
      top_feature=readLines(con, n = 1)
      top_feature=substr(top_feature, 4, nchar(top_feature))
      X_test_top_feature=data.frame(X_test[, top_feature])
      colnames(X_test_top_feature) = "value"
      p1 = plot_dist(X_test_top_feature, rownames(X_test))
      scatac_counts_plots[[fold]] = p1
      colnames(y_test) = "value"
      p2 = plot_dist(y_test, rownames(X_test))
      mut_counts_plots[[fold]] = p2
    }
    pdf("../../figures/per_fold_scatac_counts.pdf", width = 50, height = 20)
    do.call(grid.arrange, c(scatac_counts_plots, nrow=2)) 
    dev.off()
    pdf("../../figures/per_fold_mut_counts.pdf", width = 50, height = 20)
    do.call(grid.arrange, c(mut_counts_plots, nrow=2))
    dev.off()
  }
}

