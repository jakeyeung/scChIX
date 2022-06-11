# Jake Yeung
# Date of Creation: 2021-12-03
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/2-make_heatmaps_check_marker_genes.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(heatmap3)

hubprefix <- "/Users/yeung/hub_oudenaarden"
jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks

outrds.meta.colored <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/dat_meta_ordered_colorcoded.rds")
dat.ordered.lst2 <- readRDS(outrds.meta.colored)


# Load tm results ---------------------------------------------------------

inf.tm <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/tm_results_merged_lst.rds")
tm.lst <- readRDS(inf.tm)


# Get imputed -------------------------------------------------------------

dat.imputed.lst <- lapply(tm.lst, function(tm){
  t(log2(tm$topics %*% tm$terms))
})


# Load genes and  bins -------------------------------------------------------------

inf.genes <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/top_bins_genes.keeptop_250.rds")
annot.genes <- readRDS(inf.genes)

inf.annot <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/bin_annotations_from_biomart.rds"))
bins.annot <- readRDS(inf.annot)

bins2genes <- hash::hash(bins.annot$out2.df.closest$region_coord, bins.annot$out2.df.closest$gene)

inf.topicsannot <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/topics_celltypes_annotation.rds")
topicsannot <- readRDS(inf.topicsannot)

topics.hash <- hash::hash(topicsannot)
topics.hash.inv <- hash::invert(topics.hash)

# Load RNA-seq data  ------------------------------------------------------

inf.rnaseq <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/cao_rnaseq_reference_pbulk_zscore_ctypes_merged.rds")
dat.ref <- readRDS(inf.rnaseq)


# Make heatmap  -----------------------------------------------------------


clstrs.order <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "Stromal", "ConnectiveTissueProg")
names(clstrs.order) <- clstrs.order
topics.order <- sapply(clstrs.order, function(x) topics.hash.inv[[x]])

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette.sub <- cbPalette[1:length(clstrs.order)]
cols.hash <- hash::hash(clstrs.order, cbPalette.sub)

bins.vec <- unlist(lapply(annot.genes[topics.order], function(jlst) jlst$bins.keep))

bins.ctype.colvec <- lapply(clstrs.order, function(clstr){
  bin.vec <- annot.genes[[topics.order[[clstr]]]]$bins.keep
  assertthat::assert_that(length(bin.vec) >= 1)
  jcol <- cols.hash[[clstr]]
  colvec <- rep(jcol, length(bin.vec))
  names(colvec) <- bin.vec
  return(colvec)
}) %>%
  unlist() 

genes.ctype.colvec <- lapply(clstrs.order, function(clstr){
  gene.vec <- annot.genes[[topics.order[[clstr]]]]$genes.keep
  assertthat::assert_that(length(gene.vec) >= 1)
  jcol <- cols.hash[[clstr]]
  colvec <- rep(jcol, length(gene.vec))
  names(colvec) <- gene.vec
  return(colvec)
}) %>%
  unlist()
gene2col <- hash::hash(sapply(names(genes.ctype.colvec), function(x) strsplit(x, split = "\\.")[[1]][[2]]), genes.ctype.colvec)


# order cells by celltypes
dat.ordered.lst <- lapply(jmarks, function(jmark){
  dat.ordered <- dat.ordered.lst2[[jmark]] %>%
    ungroup() %>%
    filter(mark == jmark) %>%
    mutate(cluster = gsub("Precusor", "Precursor", cluster)) %>%
    mutate(cluster = factor(cluster, levels = clstrs.order)) %>%
    arrange(cluster) %>%
    rowwise() 
  dat.ordered$colorcode <- sapply(as.character(dat.ordered$cluster), function(x) cols.hash[[x]])
  return(dat.ordered)
})

# write metadata with clustercol


# add cuts.total
dat.cuts_total.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots/metadata_K36-K9me3.", jmark, ".2021-08-09.txt"))
  assertthat::assert_that(file.exists(inf.tmp))
  dat.cuts_total <- fread(inf.tmp) %>%
    dplyr::select(cell, cuts_total)
})

dat.ordered.pretty.lst <- lapply(jmarks, function(jmark){
  dat.ordered <- dat.ordered.lst[[jmark]] %>%
    left_join(., dat.cuts_total.lst[[jmark]])
  dat.ordered$clustercol <- dat.ordered$colorcode
  dat.ordered$colorcode <- NULL
  return(dat.ordered)
})
c

# write tables
outdir <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata"
for (jmark in jmarks){
  print(jmark)
  out.tmp <- file.path(outdir, paste0("metadata_cell_cluster_with_clustercol.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(dat.ordered.pretty.lst[[jmark]], file = out.tmp, quote = FALSE, sep = "\t")
}


for (jmark in jmarks){
  print(jmark)
  out.tmp <- file.path(outdir, paste0("metadata_cell_cluster_with_clustercol.singles_only.", jmark, ".", Sys.Date(), ".txt"))
  jsub <- subset(dat.ordered.pretty.lst[[jmark]], type == "single")
  fwrite(jsub, file = out.tmp, quote = FALSE, sep = "\t")
}


