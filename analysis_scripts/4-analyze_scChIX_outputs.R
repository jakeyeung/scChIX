# Jake Yeung
# Date of Creation: 2021-07-06
# File: ~/projects/scChIX/analysis_scripts/4-analyze_scChIX_outputs.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 15
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

# Load LDA  ---------------------------------------------------------------

jmark1 <- "K36"
jmark2 <- "K9m3"

jmarks <- c(jmark1, jmark2)
names(jmarks) <- jmarks

inf1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.Robj"
inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K9m3.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K9m3.K-30.Robj"
load(inf1, v=T)
out.lda.lst <- list()
out.lda.lst[[jmark1]] <- out.lda

load(inf2, v=T)
out.lda.lst[[jmark2]] <- out.lda

library(topicmodels)

tm.result.lst <- lapply(jmarks, function(jmark){
  tm.result <- posterior(out.lda.lst[[jmark]])
})

dat.umap.lst <- lapply(jmarks, function(jmark){
  jdat <- DoUmapAndLouvain(tm.result.lst[[jmark]]$topics, jsettings = jsettings)
})

dat.umap.annot.lst <- lapply(dat.umap.lst, function(jdat){
  jdat$plate <- sapply(jdat$cell, function(jcell) ClipLast(jcell, jsep = "_"))
  jdat$experi <- sapply(jdat$cell, function(jcell) ClipLast(jcell, jsep = "-"))
  return(jdat)
})

m.lst <- lapply(dat.umap.annot.lst, function(jsub){
  ggplot(jsub, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~experi) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


m.lst <- lapply(dat.umap.annot.lst, function(jsub){
  ggplot(jsub, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~plate) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


# Was there batch effect ?  -----------------------------------------------

inf.check <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000/lda_outputs.count_tables.50000.K36-K9m3.2021-06-29.K-30.binarize.FALSE/ldaOut.count_tables.50000.K36-K9m3.2021-06-29.K-30.Robj"
load(inf.check, v=T)

out.lda.check <- out.lda
tm.result.check <- posterior(out.lda.check)

dat.umap.check <- DoUmapAndLouvain(tm.result.check$topics, jsettings) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(plate, jsep = "-"))

ggplot(dat.umap.check, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# inf1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/Gastru_Unmixed_DblMark.K36.RData"
# inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/Gastru_Unmixed_DblMark.K9m3.RData"
# load(inf1, v=T)



# load(inf2, v=T)
