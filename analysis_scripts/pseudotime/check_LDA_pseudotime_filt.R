# Jake Yeung
# Date of Creation: 2021-08-06
# File: ~/projects/scChIX/analysis_scripts/pseudotime/check_LDA_pseudotime_filt.R
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
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load data  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_var_filtered_manual2nocenter.pseudotimefilt_K36_K9m3_K36-K9m3/lda_outputs.countmat_K9me3_pseudotime_filt.2021-08-05.K-30.binarize.FALSE/ldaOut.countmat_K9me3_pseudotime_filt.2021-08-05.K-30.Robj")
load(inf, v=T)


tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check pseudtime  ---------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.pseudo <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_pseudotime/metadata_pseudotime.2021-08-05.txt")
dat.pseudo <- fread(inf.pseudo)

cell2ptime.check <- hash::hash(dat.pseudo$cell, dat.pseudo$ptime)
cell2stage <- hash::hash(dat.pseudo$cell, dat.pseudo$stage)
dat.umap$ptime.check <- sapply(dat.umap$chell, function(x) cell2ptime.check[[x]])
dat.umap$stage <- sapply(dat.umap$cell, function(x) cell2stage[[x]])

ggplot(dat.umap, aes(x = umap1, y = umap2, color = ptime.check)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap, aes(x = umap1, y = umap2, color = stage)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~stage) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pseudo, aes(x = umap1, y = umap2, color = stage)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~stage) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = stage, y = umap1)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = stage, y = ptime.check)) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Refit pseudotime  -------------------------------------------------------



