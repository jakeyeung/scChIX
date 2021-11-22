# Jake Yeung
# infer_celltype_from_LDA_late_gastrulation.R
# 2021-07-17
# DESCRIPTION
# # 
# #     Load gastrulation LDA and infer cell type from topic weights
# 
# FOR HELP
# 
#     Rscript infer_celltype_from_LDA_late_gastrulation.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2021-07-17
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

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

library(assertthat)
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-ldainputvec', metavar='INFILE .RData or .rds', nargs="+",
                                            help='Input LDAs')
parser$add_argument('-marks', metavar='ANTIBODIES .RData or .rds', nargs="+",
                                            help='Input marks (order matches ldainputvec')
parser$add_argument('-pseudobulk_reference', metavar='Input reference .RData or .rds', default="/hpc/hub_oudenaarden/jyeung/data/public_data/CaoPijuana_merged_batch_cor.2019-12-03.RData",
                                            help='Input pseudobulk reference data .RData containing dat.mat.filt.batchcor or .rds')
parser$add_argument('-keeptop', metavar='TopNFeatures', default=150, type="integer",
                                            help='Keep top for N terms for comparing with reference')
parser$add_argument('-gene2tss_reference', metavar='BED', default="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed",
                                            help='Input BED file to link TSS location to gene')
parser$add_argument('-outprefix', metavar='OUTPREFIX',
                                            help='Output PREFIX. Add mark and pdf afterwards')
parser$add_argument('--calculate_var', action="store_true", default=FALSE, 
                                            help='Calculate var, default FALSE')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

keeptop <- args$keeptop
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
jinf.tss <- args$gene2tss_reference

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette <- rep(cbPalette, 10)  # repeats in case many clusters identified

inf.ref <- args$pseudobulk_reference

infs <- args$ldainputvec
assertthat::assert_that(all(file.exists(infs)))

jmarks <- args$marks
names(jmarks) <- jmarks
names(infs) <- jmarks

jstr <- paste(jmarks, collapse = "_")

# outpdf <- args$outpdf
outprefix <- args$outprefix

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load reference data to get cell types  ----------------------------------

if (endsWith(inf.ref, ".RData") | endsWith(inf.ref, ".Robj")){
  load(inf.ref, v=T)
  dat.mat.filt.batchcor <- dat.mat.filt  # no need to batch corect if we only take Shendure samples
} else if (endsWith(inf.ref, ".rds")) {
  dat.mat.filt.batchcor <- readRDS(inf.ref)
}

# keep only celltypes that start with number (shendure more late stage?)
cnames.keep <- grepl("^[[:digit:]]+", colnames(dat.mat.filt.batchcor))
dat.mat.filt.batchcor <- dat.mat.filt.batchcor[, cnames.keep]
dat.mat.filt.batchcor <- t(scale(t(dat.mat.filt.batchcor), center = TRUE, scale = TRUE))

# check batch
pca.public <- prcomp(dat.mat.filt.batchcor, center = TRUE, scale. = TRUE)
dat.pca.public <- data.frame(celltype = rownames(pca.public$x), pca.public$x, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(is.late = grepl("^[[:digit:]]+", celltype))

# m.public <- ggplot(dat.pca.public, aes(x = PC1, y = PC2, color = is.late)) +
#   geom_point() +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


genes.orig <- sapply(rownames(dat.mat.filt.batchcor), function(x) strsplit(x, "\\.")[[1]][[1]])
genes.annot <- JFuncs::EnsemblGene2Gene(gene.list = genes.orig, return.original = TRUE)
names(genes.annot) <- genes.orig

rownames(dat.mat.filt.batchcor) <- make.names(genes.annot, unique = TRUE)

dat.norm.df <- tidyr::gather(data.frame(gene = rownames(dat.mat.filt.batchcor), dat.mat.filt.batchcor), key = "celltype", value = "counts", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(counts, center = TRUE, scale = TRUE))

tm.result.lst <- lapply(infs, function(inf){
  if (endsWith(inf, ".RData") | endsWith(inf, ".Robj")){
    load(inf, v=T)  # out.lda
  } else if (endsWith(inf, ".rds")){
    out.lda <- readRDS(inf)
  } else {
    stop("Must be .RData or .rds")
  }
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(tm.result)
})



dat.umap.lst <- lapply(tm.result.lst, function(tm.result){
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)
  return(dat.umap)
})

if (args$calculate_var){
  dat.impute.log.lst <- lapply(tm.result.lst, function(tm.result){
    return(t(log2(tm.result$topics %*% tm.result$terms)))
  })
  dat.var.lst <- lapply(dat.impute.log.lst, function(dat.impute.log){
    return(CalculateVarAll(dat.impute.log, jchromos))
  })
  dat.umap.var.lst <- lapply(jmarks, function(jmark){
    print("debug 1")
    print(head(dat.umap.lst[[jmark]]))
    print(head(dat.var.lst[[jmark]]))
    dat.umap.var <- left_join(dat.umap.lst[[jmark]], dat.var.lst[[jmark]])
    return(dat.umap.var)
  })
  m.var.lst <- lapply(jmarks, function(jmark){
    print("debug 2")
    print(head(dat.umap.var.lst[[jmark]]))
    m <- ggplot(dat.umap.var.lst[[jmark]], aes(x = umap1, y = umap2, color = log(cell.var.within.sum.norm))) +
      geom_point() +
      theme_bw() +
      ggtitle(paste(jmark, "Color by Var")) +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
      scale_color_viridis_c()
    return(m)
  })
}

# Plot  -------------------------------------------------------------------


m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste(jmark, "from 50kb bins")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette)
  return(m)
})




# Get celltypes by looking at topics  -------------------------------------

# plot topics and merge with reference data

# H3K36me3 only: look at topics 50kb and assign each gene to nearest bin
# let's do TSS-TES maybe it's easier??




coords <- lapply(tm.result.lst, function(x){
  colnames(x$terms)
}) %>%
  unlist() %>%
  unique()


coords.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = coords, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

coords.annot$regions.annotated$regions_coord2 <- make.names(coords.annot$regions.annotated$region_coord)
coords.annot$out2.df$regions_coord2 <- make.names(coords.annot$out2.df$region_coord)

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
  terms.filt.tmp.merge <- left_join(terms.filt.tmp, coords.annot$out2.df, by = c("term" = "regions_coord2")) %>%
    group_by(topic) %>%
    arrange(desc(weight)) %>%
    mutate(rnk = rank(-weight))
  print(head(terms.filt.tmp.merge))

  outpdf <- paste0(outprefix, ".", jmark, ".pdf")
  pdf(outpdf, useDingbats = FALSE)

  # print(m.public)
  JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)

  if (args$calculate_var){
    # write umap with var
    print(m.var.lst)
  }

  print(m1)
  print(m2)
  print(m3)
  print(m4)
  print(m5)

  for (jtop in topics.ordered.tmp$topic){
    print(jtop)
    # m.umap <- PlotXYWithColor(dat.umap.withtopics.tmp, xvar = "umap1", yvar = "umap2", cname = jtop) + scale_color_viridis_c()
    m.umap <- ggplot(dat.umap.withtopics.tmp, aes_string(x = "umap1",  y = "umap2", color = jtop)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_viridis_c()

    top.genes <- subset(terms.filt.tmp.merge, topic == jtop & rnk <= keeptop)$gene
    assertthat::assert_that(length(top.genes) > 0)

    jsub <- subset(dat.norm.df, gene %in% top.genes)
    jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
    jlevels <- as.character(jsub.sorted.summarised$celltype)
    jsub$celltype <- factor(jsub$celltype, levels = jlevels)
    m.exprs <- ggplot(jsub,
                      aes(x = celltype , y = zscore)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
      ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
    print(m.umap)
    print(m.exprs)
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
}

