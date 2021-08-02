# Jake Yeung
# Date of Creation: 2021-07-14
# File: ~/projects/scChIX/analysis_scripts/2-check_LDA_outputs_cleaned_LessStringent2_OtherMarks.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(ggrepel)
library(scChIX)

hubprefix <- "/home/jyeung/hub_oudenaarden"
# source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/celltyping_MergedDataNoQuantNorm"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/celltyping_MergedDataNoQuantNorm_ShendureOnly_RenormScale_StricterAnnots.TopBinsFilt"
dir.create(outdir)


# Load  -------------------------------------------------------------------


jdate <- "2021-07-14"
jmarks <- c("K36", "K27", "K36-K27")

# jdate <- "2021-07-14"
# jmarks <- c("K36", "K9m3", "K36-K9m3")

names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")


dname <- "filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt_TopBinsFilt"
dname2 <- "countmat_featurefilt"
jsuffix <- paste0("50000_", dname, "_", jstr)


# jmark <- jmarks[[1]]

infs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_", jsuffix, "/lda_outputs.", dname2, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.", dname2, ".", jmark, ".", jdate, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})



tm.result.lst <- lapply(infs, function(inf){
  load(inf, v=T)  # out.lda
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(tm.result)
})
dat.umap.lst <- lapply(tm.result.lst, function(tm.result){
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)
  return(dat.umap)
})

# Plot  -------------------------------------------------------------------

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste(jmark, "from 50kb bins")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette)
  return(m)
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)


# Load reference data to get cell types  ----------------------------------

inf.ref <- file.path(hubprefix, "jyeung/data/public_data/CaoPijuana_merged_batch_cor.2019-12-03.RData")
load(inf.ref, v=T)

# dat.mat.filt.batchcor <- t(dat.mat.filt.batchcor)
dat.mat.filt.batchcor <- dat.mat.filt

# keep only celltypes that start with number (shendure more late stage?)
cnames.keep <- grepl("^[[:digit:]]+", colnames(dat.mat.filt.batchcor))
dat.mat.filt.batchcor <- dat.mat.filt.batchcor[, cnames.keep]
# mutate(is.late = grepl("^[[:digit:]]+", celltype))
# renormalize?
dat.mat.filt.batchcor <- t(scale(t(dat.mat.filt.batchcor), center = TRUE, scale = TRUE))


# check batch
pca.public <- prcomp(dat.mat.filt.batchcor, center = TRUE, scale. = TRUE)
dat.pca.public <- data.frame(celltype = rownames(pca.public$x), pca.public$x, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(is.late = grepl("^[[:digit:]]+", celltype))

ggplot(dat.pca.public, aes(x = PC1, y = PC2, color = is.late)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# # do quant norm again??? No
#
# boxplot(dat.mat.filt.batchcor)
# cnames.before <- colnames(dat.mat.filt.batchcor)
# rnames.before <- rownames(dat.mat.filt.batchcor)
# dat.mat.filt.batchcor <- preprocessCore::normalize.quantiles(dat.mat.filt.batchcor, copy=TRUE)
# colnames(dat.mat.filt.batchcor) <- cnames.before
# rownames(dat.mat.filt.batchcor) <- rnames.before

# dat.mat.filt.batchcor <- preprocessCore::normalize.quantiles(dat.mat.filt.batchcor, copy = TRUE)
# colnames(dat.mat.filt.batchcor)

genes.orig <- sapply(rownames(dat.mat.filt.batchcor), function(x) strsplit(x, "\\.")[[1]][[1]])
genes.annot <- JFuncs::EnsemblGene2Gene(gene.list = genes.orig, return.original = TRUE)
names(genes.annot) <- genes.orig

rownames(dat.mat.filt.batchcor) <- make.names(genes.annot, unique = TRUE)

# boxplot(dat.mat.filt.batchcor)

dat.norm.df <- tidyr::gather(data.frame(gene = rownames(dat.mat.filt.batchcor), dat.mat.filt.batchcor), key = "celltype", value = "counts", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(counts, center = TRUE, scale = TRUE))


# Get celltypes by looking at topics  -------------------------------------

# plot topics and merge with reference data

# H3K36me3 only: look at topics 50kb and assign each gene to nearest bin
# let's do TSS-TES maybe it's easier??


keeptop <- 150
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/celltyping"
# dir.create(outdir)


jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"

coords <- lapply(tm.result.lst, function(x){
  colnames(x$terms)
  # sapply(rownames(x$dat.raw.pbulk), function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
}) %>%
  unlist() %>%
  unique()

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


# coords.makenames <- make.names(coords)
# coords.makenames <- gsub(pattern = "\\:", "\\.", coords)
# coords.makenames <- gsub(pattern = "\\-", "\\.", coords.makenames)
coords.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = coords, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

coords.annot$regions.annotated$regions_coord2 <- make.names(coords.annot$regions.annotated$region_coord)
coords.annot$out2.df$regions_coord2 <- make.names(coords.annot$out2.df$region_coord)
# head(coords.annot$regions.annotated)

# coords.annot.lst <- lapply(coords.lst, function(coords){
# })

for (jmark in jmarks){
  print(jmark)

  topics.ordered.tmp <- OrderTopicsByEntropy(tm.result = tm.result.lst[[jmark]])



  # plot topic loadings to each UMAP
  dat.topics.tmp <- data.frame(cell = rownames(tm.result.lst[[jmark]]$topics), tm.result.lst[[jmark]]$topics, stringsAsFactors = FALSE)
  dat.umap.withtopics.tmp <- left_join(dat.umap.lst[[jmark]], dat.topics.tmp)



  # add stages
  dat.umap.withtopics.tmp$stage <- sapply(dat.umap.withtopics.tmp$cell, function(cell) StageToNumeric(GetStage(PreprocessSamp(cell))))

  # get plates
  dat.umap.withtopics.tmp$plate <- sapply(dat.umap.withtopics.tmp$cell, function(x) ClipLast(x, jsep = "_"))

  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

  m1 <- ggplot(dat.umap.withtopics.tmp, aes(x = umap1, y = umap2, color = plate)) +
    scale_color_manual(values = cbPalette) +
    geom_point() + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

  m2 <- ggplot(dat.umap.withtopics.tmp, aes(x = umap1, y = umap2, color = as.character(stage))) +
    scale_color_manual(values = cbPalette) +
    geom_point() + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

  m3 <- ggplot(dat.umap.withtopics.tmp, aes(x = umap1, y = umap2, color = plate)) +
    scale_color_manual(values = cbPalette) +
    facet_wrap(~plate) +
    geom_point() + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

  m4 <- ggplot(dat.umap.withtopics.tmp, aes(x = umap1, y = umap2, color = as.character(stage))) +
    scale_color_manual(values = cbPalette) +
    facet_wrap(~stage) +
    geom_point() + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

  m5 <- ggplot(dat.umap.withtopics.tmp, aes(x = umap1, y = umap2, color = as.character(louvain))) +
    scale_color_manual(values = cbPalette) +
    geom_point() + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")




  terms.filt.tmp <- data.frame(topic = rownames(tm.result.lst[[jmark]]$terms), as.data.frame(tm.result.lst[[jmark]]$terms)) %>%
    tidyr::gather(key = "term", value = "weight", -topic) %>%
    rowwise()
  # terms.filt.tmp.merge <- left_join(terms.filt.tmp, coords.annot$out2.df, by = c("term" = "regions_coord2"))
  terms.filt.tmp.merge <- left_join(terms.filt.tmp, coords.annot$out2.df, by = c("term" = "regions_coord2")) %>%
    # mutate(gene = ) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = rank(-weight))
  print(head(terms.filt.tmp.merge))

  outpdf <- file.path(outdir, paste0("bins_50kb_", jmark, "_celltyping_topics.", jsuffix, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)

  print(m1)
  print(m2)
  print(m3)
  print(m4)
  print(m5)

  for (jtop in topics.ordered.tmp$topic){
    print(jtop)
    # i <- strsplit(jtop, "_")[[1]][[2]]
    m.umap <- PlotXYWithColor(dat.umap.withtopics.tmp, xvar = "umap1", yvar = "umap2", cname = jtop) + scale_color_viridis_c()

    top.genes <- subset(terms.filt.tmp.merge, topic == jtop & rnk <= keeptop)$gene
    assertthat::assert_that(length(top.genes) > 0)

    jsub <- subset(dat.norm.df, gene %in% top.genes)
    jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
    jlevels <- as.character(jsub.sorted.summarised$celltype)
    jsub$celltype <- factor(jsub$celltype, levels = jlevels)
    m.exprs <- ggplot(jsub,
                      aes(x = celltype , y = zscore)) +
      geom_boxplot(outlier.shape = NA) +
      # geom_violin() +
      geom_jitter(width = 0.1, size = 0.5) +
      # geom_line() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
      ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
    print(m.umap)
    print(m.exprs)
    # plot top 150 genes?
    jsub.terms <- subset(terms.filt.tmp.merge, topic == jtop & rnk < keeptop) %>%
      ungroup() %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
    m.top <- jsub.terms %>%
      # mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
      ggplot(aes(x = term, y = log10(weight), label = gene)) +
      geom_point(size = 0.25) +
      theme_bw(8) +
      # geom_text_repel(size = keeptop / 150, segment.size = 0.1, segment.alpha = 0.25) +
      # theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = keeptop / 200)) +
      geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtop))
    print(m.top)
  }
  dev.off()
}





