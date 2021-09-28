# Jake Yeung
# Date of Creation: 2021-07-09
# File: ~/projects/scChIX/analysis_scripts/3b-clean_up_plate_effects_dbl_only.R
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

library(mixtools)

library(ggrepel)

source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")

inmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_NN_15_check_plates"

outpdf <- file.path(inmain, paste0("clean_up_clusters_check.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Load data ---------------------------------------------------------------


jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks

count.mat.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("countmat_output.", jmark, ".2021-07-07.rds")
  fpath <- file.path(inmain, fname)
  assertthat::assert_that(file.exists(fpath))
  readRDS(fpath)
})

dat.umap.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("celltyping_output.", jmark, ".2021-07-07.rds")
  fpath <- file.path(inmain, fname)
  assertthat::assert_that(file.exists(fpath))
  dat.umap.tmp <- readRDS(fpath) %>%
    rowwise() %>%
    mutate(cluster = paste("cluster", cluster, sep = ""))
})


lda.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("lda_output.", jmark, ".2021-07-07.rds")
  fpath <- file.path(inmain, fname)
  assertthat::assert_that(file.exists(fpath))
  print(fpath)
  readRDS(fpath)
})


# Clean up plate effects  -------------------------------------------------

bad.plate.dbl <- "E9p5-CB6-K36-K9m3-190409-2"

ggplot(dat.umap.lst$`K36-K9m3`, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lst$`K36-K9m3` %>% filter(plate != bad.plate.dbl), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Clean up K36 ------------------------------------------------------------

bad.clst.k36 <- "cluster13"
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.lst$`K36`, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lst$`K36` %>% filter(cluster != bad.clst.k36), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Clean up K9me3 ----------------------------------------------------------

bad.clst.k9 <- "cluster7"
ggplot(dat.umap.lst$`K9m3`, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.lst$`K9m3` %>% filter(cluster != bad.clst.k9), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(~plate) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Write output  -----------------------------------------------------------

# clean up dbl, K36, but not K9

outdir1 <- file.path(inmain, "dbl_cleaned")
dir.create(outdir1)

# write new umaps

dat.umap.filt.lst1 <- list(`K36` = dat.umap.lst$`K36`,
                          `K9m3` = dat.umap.lst$K9m3,
                          `K36-K9m3` = dat.umap.lst$`K36-K9m3` %>% filter(plate != bad.plate.dbl))



cells.keep1 <- unlist(lapply(dat.umap.filt.lst1, function(x) x$cell))
# write new count tables
count.mat.lst1 <- lapply(count.mat.lst, function(jmat){
  cols.keep.tmp <- colnames(jmat) %in% cells.keep1
  jmat[, cols.keep.tmp]
})


jdat.out1 <- lapply(jmarks, function(jmark){
  outtmp.name <- paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds")
  jdat <- dat.umap.filt.lst1[[jmark]]
  assertthat::assert_that(nrow(jdat) > 0)
  assertthat::assert_that(ncol(jdat) > 0)
  print(dim(jdat))
  saveRDS(jdat, file = file.path(outdir1, outtmp.name))
  return(jdat)
})

jmat.out1 <- lapply(jmarks, function(jmark){
  outtmp.name <- paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds")
  jmat <- count.mat.lst1[[jmark]]
  assertthat::assert_that(nrow(jmat) > 0)
  assertthat::assert_that(ncol(jmat) > 0)
  print(dim(jmat))
  saveRDS(jmat, file = file.path(outdir1, outtmp.name))
  return(jmat)
})


# write same lda
lapply(jmarks, function(jmark){
  fname <- paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds")
  fpath <- file.path(outdir1, fname)
  # assertthat::assert_that(!file.exists(fpath))
  print(fpath)
  saveRDS(lda.lst[[jmark]], file = fpath)
})


dev.off()
