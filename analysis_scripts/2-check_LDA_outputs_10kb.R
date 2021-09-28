# Jake Yeung
# Date of Creation: 2021-08-12
# File: ~/projects/scChIX/analysis_scripts/2-check_LDA_outputs_10kb.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks

infs.rdata <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_10000_var_filtered_manual2nocenternoE8unifyK36_K36_K27_K36-K27/lda_outputs.countmat_var_filt.", jmark, ".2021-08-11.K-30.binarize.FALSE/ldaOut.countmat_var_filt.", jmark, ".2021-08-11.K-30.Robj"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

out.lda.lst <- lapply(infs.rdata, function(jinf){
  load(jinf, v=T)
  return(out.lda)
})

tm.result.lst <- lapply(out.lda.lst, function(jout){
  posterior(jout)
})

dat.umap.lst <- lapply(tm.result.lst, function(tm){
  DoUmapAndLouvain(tm$topics, jsettings)
})

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

m.lst



