# Jake Yeung
# Date of Creation: 2021-10-08
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/5-recluster_cells_H3K9me3.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(topicmodels)

library(hash)
library(igraph)
library(umap)


# hubprefix <- "/Users/yeung/hub_oudenaarden"
# Recluster cells and show pseudotime ?  ----------------------------------

# load reclustering 

# inf <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/downstream/LDA_outputs/LDA_output.count_mat_gastru_H3K9me3_no_blood.2021-10-04.Robj")
inf <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/LDA_output.count_mat_gastru_H3K9me3_no_blood.2021-10-04.Robj"
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

dat.cuts.total <- data.frame(cell = colnames(count.mat), cuts_total = colSums(count.mat), stringsAsFactors = FALSE)

tm.result <- posterior(out.lda)



# Load metadata -----------------------------------------------------------

# inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/dat_meta_ordered_colorcoded.rds")
inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/metadata_cleaned.2021-10-04.rds"
dat.meta.lst <- readRDS(inf.meta)
dat.meta.long <- dat.meta.lst %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(plate = scchicFuncs::ClipLast(x = cell, jsep = "_")) %>%
  filter(mark == "K9m3") %>%
  left_join(., dat.cuts.total)

# Do umap  ----------------------------------------------------------------

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 5
jsettings[["min_dist"]] <- 0.5
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 1

nn.vec <- c(50)
mindist.vec <- c(1)

dat.umap.params <- expand.grid(nn.vec, mindist.vec) %>%
  rowwise() %>%
  mutate(nn_mindist = paste(Var1, Var2, sep = "_"))

dat.umap.params.lst <- apply(dat.umap.params, MARGIN = 1, FUN = function(jrow){
  jpair <- jrow[[3]]
  jnn <- strsplit(jpair, split = "_")[[1]][[1]]
  jmindist <- strsplit(jpair, split = "_")[[1]][[2]]
  
  jsettings <- umap.defaults
  jsettings[["n_neighbors"]] <- as.integer(jnn)
  jsettings[["min_dist"]] <- as.numeric(jmindist)
  jsettings[["random_state"]] <- 123
  jsettings[["spread"]] <- 1
  
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings) %>%
    left_join(., subset(dat.meta.long, select = c(cell, cluster, colorcode, stage))) %>%
    mutate(nn = jnn, 
           mindist = jmindist, 
           pair = jpair)
  return(dat.umap)
})

dat.umap.params.long <- dat.umap.params.lst %>%
  bind_rows() 

ggplot(dat.umap.params.long, aes(x = umap1, y = umap2, color = stage)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~pair, scales = "fixed") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Get color  --------------------------------------------------------------

dat.umap <- subset(dat.umap.params.long, pair == "50_1") %>%
  arrange(desc(stage))


m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = stage)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m)

g <- ggplot_build(m)
jcols <- unlist(unique(g$data[[1]]["colour"]))
jstages <- unique(dat.umap$stage)
stage2col <- hash::hash(jstages, jcols)

dat.umap$clustercol<- sapply(as.character(dat.umap$stage), function(x) stage2col[[x]])

m.check <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = clustercol)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

JFuncs::multiplot(m, m.check)


# Fit pseudotime  ---------------------------------------------------------


library(princurve)

# get most upper left point

# dat.umap.ptime <- dat.umap %>%
#   rowwise() %>%
#   mutate(coord = umap2)

pcurve.input <- as.matrix(subset(dat.umap, select = c(umap1, umap2)))
# pcurve.output <- princurve::principal_curve(x = pcurve.input)

pca.output <- prcomp(x = t(pcurve.input), center = FALSE, scale. = FALSE)

dat.ptime <- data.frame(cell = dat.umap$cell, dim1 = -1 * pca.output$rotation[, 1], dim2 = pca.output$rotation[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.umap) %>%
  ungroup() %>%
  mutate(pseudotime = dim1 - min(dim1))  %>%
  arrange(pseudotime) %>%
  left_join(., dat.cuts.total) %>%
  rowwise() %>%
  mutate(cluster = stage)

m.ptime <- ggplot(dat.ptime, aes(x = umap1, y = umap2, color = pseudotime)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.ptime2 <- ggplot(dat.ptime %>% arrange(desc(stage)), aes(x = umap1, y = umap2, color = stage)) + 
  geom_point(size = 3) + 
  theme_bw() + 
  scale_color_viridis_d() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Get metadata for new color ----------------------------------------------

fwrite(dat.ptime, file = "/Users/yeung/data/dblchic/gastrulation/snakemake_downstream_outputs/dat_metadata_pseudotime.H3K9me3_only.txt", sep = "\t")
# saveRDS(dat.ptime, file = "/Users/yeung/data/dblchic/gastrulation/snakemake_downstream_outputs/dat_metadata_pseudotime.H3K9me3_only.rds")


# Make plots --------------------------------------------------------------

outpdf <- paste0("/Users/yeung/data/dblchic/gastrulation/H3K9me3_pseudotime_plots/H3K9me3_pseudotime_plots.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
  print(m.ptime)
  print(m.ptime2)
dev.off()


