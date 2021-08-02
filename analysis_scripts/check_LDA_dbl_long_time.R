# Jake Yeung
# Date of Creation: 2021-07-17
# File: ~/projects/scChIX/analysis_scripts/check_LDA_dbl_long_time.R
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

# jmark <- "K9m3"
# jmark <- "K36"
# jmark <- "K27"
jmark <- "K36-K27"
# jstr <- "K36_K9m3_K36-K9m3"
jstr <- "K36_K27_K36-K27"

hubprefix <- "/home/jyeung/hub_oudenaarden"
# inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_var_filtered_manual2noblood_K36_K9m3_K36-K9m3/lda_outputs.countmat_var_filt.", jmark, ".2021-07-21.K-30.binarize.FALSE/ldaOut.countmat_var_filt.", jmark, ".2021-07-21.K-30.Robj"))
inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_var_filtered_manual2noblood_", jstr, "/lda_outputs.countmat_var_filt.", jmark, ".2021-07-21.K-30.binarize.FALSE/ldaOut.countmat_var_filt.", jmark, ".2021-07-21.K-30.Robj"))
load(inf, v=T)


tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.long <- dat.umap

# annotate
dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(plate, jsep = "-"),
         stage =  strsplit(cell, "-")[[1]][[1]],
         cluster = paste("cluster", louvain, sep = "")) %>%
  dplyr::select(-louvain)



m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  ggtitle(paste(jmark)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m1 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~plate) +
  ggtitle(paste(jmark)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~experi) +
  ggtitle(paste(jmark)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m3 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~stage) +
  ggtitle(paste( jmark)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

lapply(list(m, m1, m2, m3), print)
