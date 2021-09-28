# Jake Yeung
# Date of Creation: 2021-07-11
# File: ~/projects/scChIX/analysis_scripts/3g-load_LDA_featuresfilt_cellsfilt_output.R
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

hubprefix <- "/home/jyeung/hub_oudenaarden"

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("K36", "K9m3")
names(jmarks) <- jmarks

# Load LDA  ---------------------------------------------------------------

jdate <- "2021-07-10"
infs.lda <- lapply(jmarks, function(jmark){
  inf.lda.tmp <- file.path(hubprefix,
                           paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_topfeatures_K36_genebodies_K9m3_bins_merged/lda_outputs.countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.", jdate, ".", jmark, ".K-30.binarize.FALSE/ldaOut.countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.", jdate, ".", jmark, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf.lda.tmp))
  return(inf.lda.tmp)
})

out.objs.lst <- lapply(infs.lda, function(jinf){
  load(jinf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

out.ldas.lst <- lapply(out.objs.lst, function(jout){
  jout$out.lda
})

tm.results.lst <- lapply(out.ldas.lst, function(jlda){
  tm.result <- posterior(jlda)
  AddTopicToTmResult(tm.result)
})

count.mats.lst <- lapply(out.objs.lst, function(jout){
  jout$count.mat
})


# Load metas --------------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_LDA2/dbl_cleaned")
dats.meta <- lapply(jmarks, function(jmark){
  inf.meta.tmp <- file.path(indir.meta, paste0("celltyping_output_filt.", jmark, ".2021-07-11.rds"))
  readRDS(inf.meta.tmp) %>%
    dplyr::select(-umap1, -umap2)
})


# Make clusters -----------------------------------------------------------

dat.umaps.lst <- lapply(tm.results.lst, function(jtm){
  DoUmapAndLouvain(topics.mat = jtm$topics, jsettings = jsettings)
})

dat.umaps.merge.lst <- lapply(jmarks, function(jmark){
  left_join(dat.umaps.lst[[jmark]], dats.meta[[jmark]]) %>%
    rowwise() %>%
    mutate(stage = strsplit(cell, "-")[[1]][[1]])
})

# check louvain vscluster
m.lst <- lapply(dat.umaps.merge.lst, function(jdat){
  ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

m.clst.lst <- lapply(dat.umaps.merge.lst, function(jdat){
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

dat.umaps.merge.cleaned.lst <- lapply(dat.umaps.merge.lst, function(jdat){
  # use new louvains as cluster
  jdat$cluster <- NULL
  jdat$louvain <- paste("cluster", jdat$louvain, sep = "")
  jdat <- jdat %>%
    dplyr::rename(cluster = louvain)
})

m.cleaned.lst <- lapply(dat.umaps.merge.cleaned.lst, function(jdat){
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# Save meta, countmat, LDA objects  ---------------------------------------

outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_K36genebodies_K9m3bins_merged/dbl_cleaned")
dir.create(outdir)

lapply(jmarks, function(jmark){
  metaname <- paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds")
  metapath <- file.path(outdir, metaname)
  jtmp <- dat.umaps.merge.cleaned.lst[[jmark]]
  print(dim(jtmp))
  saveRDS(jtmp, file = metapath)
})

lapply(jmarks, function(jmark){
  metaname <- paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds")
  metapath <- file.path(outdir, metaname)
  jtmp <- count.mats.lst[[jmark]]
  print(dim(jtmp))
  saveRDS(jtmp, file = metapath)
})

lapply(jmarks, function(jmark){
  metaname <- paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds")
  metapath <- file.path(outdir, metaname)
  jtmp <- out.ldas.lst[[jmark]]
  print(dim(jtmp))
  saveRDS(jtmp, file = metapath)
})

# write pdfs
lapply(jmarks, function(jmark){
  outpdfname <- paste0("plots.", jmark, ".", Sys.Date(), ".pdf")
  outpdfpath <- file.path(outdir, outpdfname)
  pdf(file = outpdfpath, useDingbats = FALSE)
  print(m.cleaned.lst[[jmark]])
  dev.off()
})


# write meta and countmat for dbl

# load dbl countmat
jmarkdbl <- "K36-K9m3"
inf.count.mat.dbl <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/coords_filtered_K36_K9m3_merged_with_dbl/countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.2021-07-11.", jmarkdbl, ".rds"))
assertthat::assert_that(file.exists(inf.count.mat.dbl))
count.mat.dbl <- readRDS(inf.count.mat.dbl)

# load dbl meta
inf.meta.dbl <- file.path(indir.meta, paste0("celltyping_output_filt.", jmarkdbl, ".2021-07-11.rds"))
assertthat::assert_that(file.exists(inf.meta.dbl))
dat.meta.dbl <- readRDS(inf.meta.dbl)

metapath.dbl <- file.path(outdir, paste0("celltyping_output_filt.", jmarkdbl, ".", Sys.Date(), ".rds"))
countpath.dbl <- file.path(outdir, paste0("countmat_output_filt.", jmarkdbl, ".", Sys.Date(), ".rds"))
ldapath.dbl <- file.path(outdir, paste0("lda_output_filt.", jmarkdbl, ".", Sys.Date(), ".rds"))

inf.lda.dbl <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_topfeatures_K36_genebodies_K9m3_bins_merged_with_dbl/lda_outputs.countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.2021-07-11.", jmarkdbl, ".K-30.binarize.FALSE/ldaOut.countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.2021-07-11.", jmarkdbl, ".K-30.Robj"))
load(inf.lda.dbl, v=T)

out.lda.dbl <- out.lda

saveRDS(count.mat.dbl, file = countpath.dbl)
saveRDS(dat.meta.dbl, file = metapath.dbl)
saveRDS(out.lda.dbl, file = ldapath.dbl)

