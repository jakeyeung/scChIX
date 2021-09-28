# Jake Yeung
# Date of Creation: 2021-07-12
# File: ~/projects/scChIX/analysis_scripts/check_K36_K9me3_dbl_cleaned_clusters.R
# Why K9me3 cluster6 has no unmixed cells?

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(scChIX)
library(irlba)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load meta ---------------------------------------------------------------


dats.meta <- lapply(jmarks, function(jmark){
  # inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_cleaned_K36genebodies_K9m3bins_merged/dbl_cleaned/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.Robj")
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_K36genebodies_K9m3bins_merged/dbl_cleaned/celltyping_output_filt.", jmark, ".2021-07-11.rds"))
  assertthat::assert_that(file.exists(inf.meta))
  readRDS(inf.meta)
})

lapply(dats.meta, dim)

ggplot(dats.meta$K9m3 %>% filter(stage == "E9p5"), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~plate) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load raw mat  -----------------------------------------------------------


countmats <- lapply(jmarks, function(jmark){
  inf.countmat <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_K36genebodies_K9m3bins_merged/dbl_cleaned/countmat_output_filt.", jmark, ".2021-07-11.rds"))
  assertthat::assert_that(file.exists(inf.countmat))
  readRDS(inf.countmat)
})

countmats.merge <- do.call(cbind, countmats)


# Do LSI  -----------------------------------------------------------------

dat.lsi <- RunLSI(count.mat = as.matrix(countmats.merge))

dat.umap.lsi <- DoUmapAndLouvain(dat.lsi$u, jsettings = jsettings)

dat.umap.lsi.annot <- dat.umap.lsi %>%
  rowwise() %>%
  mutate(mark = strsplit(cell, split = "-")[[1]][[3]])

ggplot(dat.umap.lsi.annot, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  facet_wrap(~mark) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load raw mats before filtering ------------------------------------------

# inf.mat.orig <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_50000/K36/cbind_out/K36.countTable.binsize_50000.rds"

mats.orig <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.mat.orig <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_topfeatures_K36_genebodies_K9m3_bins/", jmark, "/cbind_out/", jmark, ".counts_tables_topfeatures_K36_genebodies_K9m3_bins.rds")
  assertthat::assert_that(file.exists(inf.mat.orig))
  mat.orig <- readRDS(inf.mat.orig)
})

mats.orig$K36[1:5, 1:5]
mats.orig$K9m3[1:5, 1:5]

dat.lsi.orig.lst <- lapply(jmarks, function(jmark){
  dat.lsi.orig <- RunLSI(count.mat = as.matrix(mats.orig[[jmark]]))
  dat.umap.lsi.orig <- DoUmapAndLouvain(dat.lsi.orig$u, jsettings = jsettings)

  dat.umap.lsi.orig.annot <- dat.umap.lsi.orig %>%
    rowwise() %>%
    mutate(stage = strsplit(cell, split = "-")[[1]][[1]])

  # label good or bad
  cells.good <- dats.meta[[jmark]]$cell

  dat.umap.lsi.orig.annot <- dat.umap.lsi.orig.annot %>%
    rowwise() %>%
    mutate(is.good = cell %in% cells.good)
})



ggplot(dat.lsi.orig.lst[[1]], aes(x = umap1, y = umap2, color = stage)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~stage) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.lsi.orig.lst[[2]], aes(x = umap1, y = umap2, color = stage)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.lsi.orig.lst[[1]], aes(x = umap1, y = umap2, color = is.good)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.lsi.orig.lst[[2]], aes(x = umap1, y = umap2, color = is.good)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# check K9me3 plate effects?

jcheck <- dat.lsi.orig.lst$K9m3 %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         experi = ClipLast(plate, jsep = "-"))


ggplot(jcheck, aes(x = umap1, y = umap2, color = stage)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~plate) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jcells.island <- subset(dats.meta$K9m3, umap2 > 3)$cell

ggplot(jcheck %>% mutate(is.blood = cell %in% jcells.island), aes(x = umap1, y = umap2, color = is.blood)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~plate) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




ggplot(dats.meta$K9m3, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~plate) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
