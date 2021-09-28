# Jake Yeung
# Date of Creation: 2021-09-10
# File: ~/projects/scChIX/analysis_scripts/unfixed/unfixed_10kb_LDA_downstream.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load data  --------------------------------------------------------------

jmarks <- c("K4m1", "K27m3", "K4m1_K27m3"); names(jmarks) <- jmarks

jsuffix <- "10kb_genomewide"
infs.lda <- lapply(jmarks, function(jmark){
  print(jmark)
  # inf.tmp <- paste0(hubprefix, "jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_10kb_genomewide_filt/lda_outputs.count_mat.10kb_genomewide_filt.", jmark, ".K-30.binarize.FALSE/ldaOut.count_mat.10kb_genomewide_filt.", jmark, ".K-30.Robj")
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_", jsuffix, "/lda_outputs.count_mat.", jsuffix, ".", jmark, ".K-30.binarize.FALSE/ldaOut.count_mat.", jsuffix, ".", jmark, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

out.lda.lst <- lapply(infs.lda, function(jinf){
  load(jinf, v=T)
  return(out.lda)
})

tm.lst <- lapply(jmarks, function(jmark){
  topicmodels::posterior(out.lda.lst[[jmark]])
})


# Load metas --------------------------------------------------------------

inf.rdata <- file.path(hubprefix, "jyeung/data/dblchic/from_rstudio/primetime/unfixed_louvain2/BM_UnfixedLouvain2.FinalCellClusterTable.2020-03-21.RData")
load(inf.rdata, v=T)


# Do umaps  ---------------------------------------------------------------


dat.umaps.lst <- lapply(jmarks, function(jmark){
  topics.mat <- tm.lst[[jmark]]$topics
  dat.umap.tmp <- scchicFuncs::DoUmapAndLouvain(topics.mat, jsettings)
  if (jmark != jmarks[[3]]){
    dat.umap.tmp <- dat.umap.tmp %>%
      left_join(., subset(dat.final.annots, mark == jmark, select = c(cell, cluster, plate, mark)))
  } else {
    dat.umap.tmp <- dat.umap.tmp %>%
      left_join(., subset(dat.final.annots, mark == jmarks[[1]], select = c(cell, cluster, plate, mark)))
  }
})

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umaps.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster))  +
    geom_point() +
    scale_color_manual(values = cbPalette) +
    theme_bw() +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

print(m.lst)




