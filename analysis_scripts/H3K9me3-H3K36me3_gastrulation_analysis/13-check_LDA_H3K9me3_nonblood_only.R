# Jake Yeung
# Date of Creation: 2021-12-28
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/13-check_LDA_H3K9me3_nonblood_only.R
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

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123


inf <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/lda_objs/ldaAnalysis_50000_var_filtered_manual2nocenter.pseudotimefilt.chr9binremoved_K36_K9m3_K36-K9m3/lda_outputs.countmat_K9me3_pseudotime_filt.badbin_filt.2021-12-28.K-30.binarize.FALSE/ldaOut.countmat_K9me3_pseudotime_filt.badbin_filt.2021-12-28.K-30.Robj"
load(inf, v=T)


tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

# add pseudotime
dat.umap$stage <- sapply(as.character(dat.umap$cell), function(x) strsplit(x, "-")[[1]][[1]])

dat.umap$stage <- factor(dat.umap$stage, levels = c("E9p5", "E10p5", "E11p5"))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = stage))  + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.umap, aes(x = umap2, fill = stage))  + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = stage, y = umap2))  + 
  geom_boxplot() + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





