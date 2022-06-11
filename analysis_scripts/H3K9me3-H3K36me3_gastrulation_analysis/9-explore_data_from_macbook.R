# Jake Yeung
# Date of Creation: 2021-11-16
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/9-explore_data_from_macbook.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata/metadata_flipped.rds"
dat.meta <- readRDS(inf)

inf.dbl <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/H3K36me3_H3K9me3_downstream_objs/dbl_cell_assignment_to_cluster.2021-10-22.rds"
dat.dbl <- readRDS(inf.dbl)


# Load obj do LDA original ------------------------------------------------

inf1 <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/objs_from_LDA/celltyping_output_filt.K36.rds"
inf2 <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/objs_from_LDA/celltyping_output_filt.K9m3.rds"

dat1 <- readRDS(inf1)

clst.merge <- c("cluster1", "cluster5", "cluster3")
dat2 <- readRDS(inf2) %>%
  mutate(cluster2 = cluster,
         cluster2 = ifelse(cluster %in% clst.merge, "cluster1", cluster),
         cluster2 = ifelse(cluster2 == "cluster4", "cluster3", cluster2))

cbPalette <- c("#696969",  "#f80597", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187")

ggplot(dat1, aes(x = -1 * umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  xlab("umap1") + ylab("umap2") + 
  ggtitle("H3K36me3 active mark") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat2, aes(x = umap1, y = -1 * umap2, color = cluster2)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  xlab("umap1") + ylab("umap2") + 
  ggtitle("H3K9me3 heterochromatin mark") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



