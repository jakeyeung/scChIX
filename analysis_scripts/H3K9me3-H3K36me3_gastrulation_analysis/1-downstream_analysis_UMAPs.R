# Jake Yeung
# Date of Creation: 2021-09-27
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/1-downstream_analysis_UMAPs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(igraph)
library(umap)
library(scChIX)
library(topicmodels)

library(scchicFuncs)

library(ggrepel)

hubprefix <- "/Users/yeung/hub_oudenaarden"
jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks


# Functions ---------------------------------------------------------------


MeanAcrossClusters2 <- function(count.mat, cnames.keep.lst, jfunc = rowMeans){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    if (length(cnames.keep.i) == 1){
      return(count.mat[, cnames.keep.i])
    } else {
      jfunc(count.mat[, cnames.keep.i])
    }
  })
  return(count.vecs)
}


# Load data  --------------------------------------------------------------

infs.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/projection_output.", jmark, ".RData"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

# infs.lst <- lapply(jmarks, function(jmark){
#   inf.tmp <- paste0("/Users/yeung/Dropbox/from_cluster/dblchic/gastrulation_objs/K36_K9m3_K36-K9m3/projection_output.", jmark, ".RData")
#   assertthat::assert_that(file.exists(inf.tmp))
#   return(inf.tmp)
# })


outs.lst <- lapply(infs.lst, function(jinf){
  load(jinf, v=T)
  return(list(out.lda.predict = out.lda.predict, out.objs = out.objs))
})



# Load metas  -------------------------------------------------------------

inf.meta.all <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_demux_cleaned_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/demux_cleaned_filtered_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3.2021-08-24.filt2.spread_7.single_and_dbl.txt")
# inf.meta.all <- paste0("/Users/yeung/Dropbox/from_cluster/dblchic/gastrulation_objs/K36_K9m3_K36-K9m3/metadata/from_demux_cleaned_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/demux_cleaned_filtered_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3.2021-08-24.filt2.spread_7.single_and_dbl.txt")
dat.meta.all <- fread(inf.meta.all)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ggplot(dat.meta.all, aes(x = umap1.shift, y = umap2.scale, color = cluster, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Flip K36 ----------------------------------------------------------------

dat.meta.all.flipped <- dat.meta.all %>%
  group_by(mark) %>%
  mutate(umap1.flip = ifelse(mark == "K36", umap1.scale * -1, umap1.scale), 
         umap1.shift2 = ifelse(mark == "K36", umap1.flip - 5, umap1.shift), 
         umap2.flip = ifelse(mark == "K36", umap2.scale, umap2.scale * -1)) %>%
  ungroup() %>%
  mutate(stage = factor(x = stage, levels = c("E9p5", "E10p5", "E11p5"))) %>%
  arrange(desc(stage))

ggplot(dat.meta.all.flipped, aes(x = umap1.shift2, y = umap2.flip, color = cluster, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Show pseudotime  --------------------------------------------------------

clsts.remove <- c("Erythroid", "WhiteBloodCells")

jsub <- dat.meta.all.flipped %>% filter(mark == "K9m3" & !cluster %in% clsts.remove)

ggplot(jsub,
       aes(x = umap1.shift2, y = umap2.flip, color = stage)) + 
  geom_point() + 
  scale_color_brewer(palette = "Reds") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Check cell types across stages ------------------------------------------

# jsub.sum <- jsub %>%
#   ungroup() %>%
#   mutate(umap1.round = round(umap1.shift2, 1)) %>%
#   group_by(stage, cluster) %>%
#   summarise(ncells = length(cell)) %>%
#   group_by(stage) %>%
#   mutate(nfrac = ncells / sum(ncells))

jsub.sum <- jsub %>%
  ungroup() %>%
  mutate(umap1.round = round(umap1.shift2 / 0.5, digits = 0) * 0.5) %>%
  group_by(umap1.round, cluster) %>%
  summarise(ncells = length(cell)) %>%
  group_by(umap1.round) %>%
  mutate(nfrac = ncells / sum(ncells))
print(unique(jsub.sum$umap1.round))

ggplot(jsub.sum, aes(x = umap1.round, y = nfrac, color = cluster, group = cluster)) + 
  geom_line(size = 3) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Do louvain clustering on H3K9me3 and show it does not work  -------------

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 5

# cluster of K36me3
jmark <- "K36"

wrangled.out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  tm.orig <- posterior(outs.lst[[jmark]]$out.objs$out.lda)
  tm.orig <- AddTopicToTmResult(tm.orig)
  
  tm.proj <- outs.lst[[jmark]]$out.lda.predict
  tm.proj <- AddTopicToTmResult(tm.proj)
  
  tm.topics.merged <- rbind(tm.orig$topics, tm.proj$topics)
  
  umap.out <- umap(tm.orig$topics, config = jsettings)
  dat.umap.orig <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  dat.umap.orig <- DoLouvain(topics.mat = tm.orig$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.orig)
  # add projections
  umap.out.pred.layout <- predict(umap.out, data = tm.proj$topics)
  dat.umap.pred <- data.frame(cell = rownames(umap.out.pred.layout), 
                                    umap1 = umap.out.pred.layout[, 1], 
                                    umap2 = umap.out.pred.layout[, 2], 
                                    stringsAsFactors = FALSE)
  dat.umap.pred <- DoLouvain(topics.mat = tm.proj$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.pred)
  
  dat.umap.merged <- rbind(dat.umap.orig, dat.umap.pred)
  
  # dat.umap.merged <- DoUmapAndLouvain(tm.topics.merged, jsettings)
  dat.umap.merged.annot <- left_join(dat.umap.merged, subset(dat.meta.all, mark = jmark, select = c(cell, cluster, type, stage, mark)))
  return(list(dat.umap.merged.annot = dat.umap.merged.annot, tm.topics.merged = tm.topics.merged, tm.terms = tm.orig$terms))
})

dat.umap.merged.annot.lst <- lapply(jmarks, function(jmark){
  wrangled.out.lst[[jmark]]$dat.umap.merged.annot
})

tm.topics.merged.lst <- lapply(jmarks, function(jmark){
  wrangled.out.lst[[jmark]]$tm.topics.merged
})

tm.terms.merged.lst <- lapply(jmarks, function(jmark){
  wrangled.out.lst[[jmark]]$tm.terms
})

m.louv.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.merged.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

m.clst.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.merged.annot.lst[[jmark]] %>% rowwise(), aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    # facet_wrap(~type) + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

JFuncs::multiplot(m.louv.lst$K36, m.louv.lst$K9m3, cols = 2)
JFuncs::multiplot(m.clst.lst$K36, m.clst.lst$K9m3, cols = 2)

m.clst.check.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.merged.annot.lst[[jmark]] %>% rowwise() %>% mutate(cluster = cluster == "Neurons"), aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    theme_bw() + 
    # facet_wrap(~type) + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})
print(m.clst.check.lst$K36)

# Assign clusters using topics --------------------------------------------
# K36

# 1) Erythroids (topic10)
# 2) Neuron progenitors (topic11)
# 3) Epithelial cells? (topic21)
# 4) Endothelial cells (topic17)
# 5) Cardiac muscle progenitors? (topic5)
# 6) Neural progenitors, neural tube (topic1)
# 7) Chondrocyte progenitors (topic6)
# 8) White blood cells (topic13)

topics.keep.lst <- list("topic27" = "Endothelial", 
                        "topic24" = "Erythroid",
                        "topic29" = "Epithelial", 
                        "topic25" = "Neurons", 
                        "topic26" = "ConnectiveTissueProg", 
                        "topic20" = "NeuralTubeNeuralProgs", 
                        "topic15" = "WhiteBloodCells", 
                        "topic13" = "SchwannCellPrecursor", 
                        "topic17" = "Stromal")

rename.list <- list("Stromal" = "MesenchymalProgenitors")

# topics.keep <- c("topic10", "topic11", "topic21", "topic17", "topic5", "topic1", "topic6", "topic13")
# topics.keep.name <- c("Erythroid", "NeuronProgs", "Epithelial", "Endothelial", "CardiacMuscProg", "NeuralTube", "ChondrocyteProgs", "WhiteBlood")
# names(topics.keep.name) <- topics.keep
# assertthat::assert_that(length(topics.keep) == length(topics.keep.name))

# check
# (jtop <- topics.keep[[1]])

jmarktmp <- "K36"
tm.result.merged <- list(topics = tm.topics.merged.lst[[jmarktmp]], terms = tm.terms.merged.lst[[jmarktmp]])
dat.topics.ordered <- OrderTopicsByEntropy(tm.result = tm.result.merged)

# jtop <- "topic1"
outpdf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/celltyping_check_from_merged/celltyping_", jmarktmp, ".", Sys.Date(), ".pdf"))
pdf(file = outpdf, useDingbats = FALSE)
print(m.clst.lst$K36)
for (jtop in dat.topics.ordered$topic){
  print(jtop)
  dat.topics.tmp <- data.frame(cell = rownames(tm.topics.merged.lst[[jmarktmp]]), topic.weight = tm.topics.merged.lst[[jmarktmp]][, jtop])
  dat.check <- dat.umap.merged.annot.lst[[jmarktmp]] %>%
    left_join(., dat.topics.tmp)
  
  m.check <- ggplot(dat.check, aes(x = umap1, y = umap2, color = topic.weight)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtop) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.check)
}
dev.off()


# Load reference data to get cell types  ----------------------------------

inf.ref <- file.path(hubprefix, "jyeung/data/public_data/CaoPijuana_merged_batch_cor.2019-12-03.RData")
load(inf.ref, v=T)

dat.mat.filt.batchcor <- dat.mat.filt
# keep only celltypes that start with number (shendure more late stage?)
cnames.keep <- grepl("^[[:digit:]]+", colnames(dat.mat.filt.batchcor))
dat.mat.filt.batchcor <- dat.mat.filt.batchcor[, cnames.keep]
dat.mat.filt.batchcor <- t(scale(t(dat.mat.filt.batchcor), center = TRUE, scale = TRUE))

# check batch
pca.public <- prcomp(dat.mat.filt.batchcor, center = TRUE, scale. = TRUE)
dat.pca.public <- data.frame(celltype = rownames(pca.public$rotation), pca.public$rotation, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(is.late = grepl("^[[:digit:]]+", celltype))

ggplot(dat.pca.public, aes(x = PC1, y = PC2, label = celltype)) +
  geom_point() +
  geom_text() + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


genes.orig <- sapply(rownames(dat.mat.filt.batchcor), function(x) strsplit(x, "\\.")[[1]][[1]])

# saveRDS(genes.orig, "/Users/yeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/check_genes_orig.rds")
# run on server

# library(JFuncs)
# genes.annot <- JFuncs::EnsemblGene2Gene(gene.list = genes.orig, return.original = TRUE)
# names(genes.annot) <- genes.orig
outf.genes.annot <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/genes_annotated.rds")
genes.annot <- readRDS(file = outf.genes.annot)

rownames(dat.mat.filt.batchcor) <- make.names(genes.annot, unique = TRUE)

dat.norm.df <- tidyr::gather(data.frame(gene = rownames(dat.mat.filt.batchcor), dat.mat.filt.batchcor), key = "celltype", value = "counts", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(counts, center = TRUE, scale = TRUE))

# Find genes associated with cell types  ----------------------------------

# annotate bins
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

bins.all <- colnames(tm.result.merged$terms)
jinf.tss <- "/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.annot <- AnnotateCoordsFromList(coords.vec = bins.all, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

bins.annot$out2.df.closest$regions_coord2 <- make.names(bins.annot$out2.df.closest$region_coord)

jtops <- names(topics.keep.lst); names(jtops) <- jtops

keeptop <- 150

terms.filt.tmp <- data.frame(topic = rownames(tm.result.merged$terms), as.data.frame(tm.result.merged$terms)) %>%
  tidyr::gather(key = "term", value = "weight", -topic) %>%
  rowwise()
terms.filt.tmp.merge <- left_join(terms.filt.tmp, bins.annot$out2.df, by = c("term" = "regions_coord2")) %>%
  group_by(topic) %>%
  arrange(desc(weight)) %>%
  mutate(rnk = rank(-weight))


outpdf2 <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/celltyping_check_from_merged/celltyping_check_RNAseq_show_top_terms.", jmarktmp, ".", Sys.Date(), ".pdf"))
pdf(file = outpdf2, useDingbats = FALSE)
jsub.terms.lst <- lapply(jtops, function(jtop){
  print(jtop)
  bins.vec.tmp <- sort(tm.result.merged$terms[jtop, ], decreasing = TRUE)
  bins.keep <- names(bins.vec.tmp[1:keeptop])
  top.genes <- subset(bins.annot$out2.df.closest, region_coord %in% bins.keep)$gene
  
  assertthat::assert_that(length(top.genes) > 0)
  
  jsub <- subset(dat.norm.df, gene %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
  jlevels <- as.character(jsub.sorted.summarised$celltype)
  jsub$celltype <- factor(jsub$celltype, levels = jlevels)
  # print(head(jsub))
  assertthat::assert_that(nrow(jsub) > 0)
  m.exprs <- ggplot(jsub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
    ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  print(m.exprs)
  
  # plot top 150 genes?
  jsub.terms <- subset(terms.filt.tmp.merge, topic == jtop & rnk < keeptop) %>%
    ungroup() %>%
    mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  m.top <- jsub.terms %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(8) +
    geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(paste("Top peak weights for:", jtop))
  print(m.top)
  # return topgenes
  return(jsub.terms)
})
dev.off()

# check again

keeptop <- 250
outpdf3 <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/celltyping_check_from_merged/celltyping_check_final.keeptop_", keeptop, ".", jmarktmp, ".", Sys.Date(), ".pdf"))
pdf(file = outpdf3, useDingbats = FALSE)
for (jtop in dat.topics.ordered$topic){
  print(jtop)
  dat.topics.tmp <- data.frame(cell = rownames(tm.topics.merged.lst[[jmarktmp]]), topic.weight = tm.topics.merged.lst[[jmarktmp]][, jtop])
  dat.check <- dat.umap.merged.annot.lst[[jmarktmp]] %>%
    left_join(., dat.topics.tmp)
  
  m.check <- ggplot(dat.check, aes(x = umap1, y = umap2, color = topic.weight)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtop) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.check)
  
  bins.vec.tmp <- sort(tm.result.merged$terms[jtop, ], decreasing = TRUE)
  bins.keep <- names(bins.vec.tmp[1:keeptop])
  top.genes <- subset(bins.annot$out2.df.closest, region_coord %in% bins.keep)$gene
  
  assertthat::assert_that(length(top.genes) > 0)
  
  jsub <- subset(dat.norm.df, gene %in% top.genes)
  jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
  jlevels <- as.character(jsub.sorted.summarised$celltype)
  jsub$celltype <- factor(jsub$celltype, levels = jlevels)
  # print(head(jsub))
  assertthat::assert_that(nrow(jsub) > 0)
  m.exprs <- ggplot(jsub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
    ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  print(m.exprs)
  
  # plot top 150 genes?
  jsub.terms <- subset(terms.filt.tmp.merge, topic == jtop & rnk < keeptop) %>%
    ungroup() %>%
    mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  m.top <- jsub.terms %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(8) +
    geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(paste("Top peak weights for:", jtop))
  print(m.top)
}
dev.off()


# Get gene sets from topics  ----------------------------------------------

topics.filt <- names(topics.keep.lst); names(topics.filt) <- topics.filt

top.bins.genes <- lapply(topics.filt, function(jtop){
  print(jtop)
  bins.vec.tmp <- sort(tm.result.merged$terms[jtop, ], decreasing = TRUE)
  bins.keep <- names(bins.vec.tmp[1:keeptop])
  top.genes <- subset(bins.annot$out2.df.closest, region_coord %in% bins.keep)$gene
  assertthat::assert_that(length(top.genes) > 0)
  return(list(bins.keep = bins.keep, genes.keep = top.genes))
})

# topics to celltype


# Save outputs ------------------------------------------------------------

outf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/top_bins_genes.keeptop_", keeptop, ".rds"))
saveRDS(object = top.bins.genes, file = outf)

outf.annot <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/bin_annotations_from_biomart.rds"))
saveRDS(object = bins.annot, file = outf.annot)

outf.topicsannot <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/topics_celltypes_annotation.rds"))
saveRDS(object = topics.keep.lst, outf.topicsannot)

outf.rnaseq <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/cao_rnaseq_reference_pbulk_zscore.rds"))
saveRDS(dat.norm.df, file = outf.rnaseq)

# save tm results merged
outf.tm <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/tm_results_merged_lst.rds"))
tm.result.merged.lst <- lapply(jmarks, function(jmarktmp){
  print(jmarktmp)
  tm.result.tmp <- list(topics = tm.topics.merged.lst[[jmarktmp]], terms = tm.terms.merged.lst[[jmarktmp]])
  print(head(tm.result.tmp$topics[1:5, 1:5]))
  return(tm.result.tmp)
})
saveRDS(tm.result.merged.lst, file = outf.tm)

# write meta

outf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/metadata_flipped.rds"))
saveRDS(dat.meta.all.flipped, file = outf.meta)


# Merge some celltypes together -------------------------------------------

dat.mat.filt.cao <- dat.mat.filt[, cnames.keep]

print(sort(colnames(dat.mat.filt.cao), decreasing = TRUE))

cnames <- colnames(dat.mat.filt.cao)
# merge erythroids
eryths <- cnames[grepl("eryth", cnames)]
neurons <- cnames[grepl("neur|Neural.progenitor", cnames)]
mesench <- cnames[grepl("Chond|mesench|tooth|Connective|Osteoblasts|Meso|Stromal", cnames)]
neurtube <- cnames[grepl("Neural.Tube|Radial|ligodendrocyte|Notochord|Isthmic", cnames)]

jnames <- cnames[!cnames %in% c(eryths, neurons, mesench, neurtube)]
jnames.toadd <- c("Erythroid", "NeuronalProgenitors", "MesenchymalProgenitors", "NeuralTubeAndProgs")
jnames.appended <- c(jnames, jnames.toadd)
names(jnames.appended) <- jnames.appended

cnames.toadd <- list(eryths, neurons, mesench, neurtube); names(cnames.toadd) <- jnames.toadd

cnames.keep.lst <- lapply(jnames.appended, function(jname){
  if (grepl("^[[:digit:]]+", jname)){
    return(jname)
  } else {
    return(cnames.toadd[[jname]])
  }
})

pbulk.mat <- as.data.frame(MeanAcrossClusters2(count.mat = as.matrix(dat.mat.filt.cao), cnames.keep.lst = cnames.keep.lst, jfunc = rowMeans))
pbulk.mat.scaled <- t(scale(t(pbulk.mat), center = FALSE, scale = FALSE))

genes.hash <- hash::hash(names(genes.annot), genes.annot)
rownames(pbulk.mat.scaled) <- sapply(rownames(pbulk.mat.scaled), function(x) AssignHash(x = x, jhash = genes.hash, null.fill = NA))

dat.norm.merged.df <- tidyr::gather(data.frame(gene = rownames(pbulk.mat.scaled), pbulk.mat.scaled), key = "celltype", value = "counts", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(counts, center = TRUE, scale = TRUE))

outf.rnaseq.merged <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/cao_rnaseq_reference_pbulk_zscore_ctypes_merged.rds"))
saveRDS(object = dat.norm.merged.df, file = outf.rnaseq.merged)

# Redo  -------------------------------------------------------------------

outpdf.merged <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/celltyping_check_from_merged/celltyping_check_final.ctypes_merged.keeptop_", keeptop, ".", jmarktmp, ".", Sys.Date(), ".pdf"))
pdf(file = outpdf.merged, useDingbats = FALSE)
for (jtop in dat.topics.ordered$topic){
  print(jtop)
  dat.topics.tmp <- data.frame(cell = rownames(tm.topics.merged.lst[[jmarktmp]]), topic.weight = tm.topics.merged.lst[[jmarktmp]][, jtop])
  dat.check <- dat.umap.merged.annot.lst[[jmarktmp]] %>%
    left_join(., dat.topics.tmp)
  
  m.check <- ggplot(dat.check, aes(x = umap1, y = umap2, color = topic.weight)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jtop) + 
    scale_color_viridis_c() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.check)
  
  bins.vec.tmp <- sort(tm.result.merged$terms[jtop, ], decreasing = TRUE)
  bins.keep <- names(bins.vec.tmp[1:keeptop])
  top.genes <- subset(bins.annot$out2.df.closest, region_coord %in% bins.keep)$gene
  
  assertthat::assert_that(length(top.genes) > 0)
  
  jsub <- subset(dat.norm.merged.df, gene %in% top.genes)
  jsub.sorted.summarised <- jsub %>% 
    group_by(celltype) %>% 
    summarise(zscore = median(zscore)) %>% 
    arrange(desc(zscore)) %>% 
    dplyr::select(celltype)
  jlevels <- as.character(jsub.sorted.summarised$celltype)
  jsub$celltype <- factor(jsub$celltype, levels = jlevels)
  # print(head(jsub))
  assertthat::assert_that(nrow(jsub) > 0)
  m.exprs <- ggplot(jsub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_violin() +
    geom_jitter(width = 0.1, size = 0.5) +
    # geom_line() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
    ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  print(m.exprs)
  
  # plot top 150 genes?
  jsub.terms <- subset(terms.filt.tmp.merge, topic == jtop & rnk < keeptop) %>%
    ungroup() %>%
    mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  m.top <- jsub.terms %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(8) +
    geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(paste("Top peak weights for:", jtop))
  print(m.top)
}
dev.off()



