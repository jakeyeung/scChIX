# Jake Yeung
# Date of Creation: 2021-07-14
# File: ~/projects/scChIX/analysis_scripts/3e-get_top_bins_and_genes.LessStringent2.R
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

topbins.keep <- 250
GetTopFeatureNames <- function(jrow, cnames, jtopbins.keep){
  names(jrow) <- cnames
  jrow.rank <- rank(-1 * jrow)
  jrow.rank.keep <- jrow.rank < topbins.keep
  features.keep <- names(jrow.rank.keep)[jrow.rank.keep]
  return(features.keep)
}

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load bins and TSSs ------------------------------------------------------

# jmarks.single <- c("K36", "K9m3")
# jmarks.dbl <- c("K36-K9m3")

jmarks.single <- c("K36", "K27")
jmarks.dbl <- c("K36-K27")

jmarks <- c(jmarks.single, jmarks.dbl)
names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# load the two ldas

dname <- "filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt"
inbase <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2")
inmain <- file.path(inbase, jmarks[[2]], dname)
outmain <- file.path(inbase, "coords_filtered")
dir.create(outmain)
outdir <- file.path(outmain, jmarks[[2]])
dir.create(outdir)


infs.lda <- lapply(jmarks, function(jmark){
  inf.lda.tmp <- file.path(inmain, paste0("lda_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.lda.tmp))
  return(inf.lda.tmp)
})

infs.countmat <- lapply(jmarks, function(jmark){
  inf.countmat.tmp <- file.path(inmain, paste0("countmat_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.countmat.tmp))
  return(inf.countmat.tmp)
})


lda.objs <- lapply(infs.lda, readRDS)

count.mats <- lapply(infs.countmat, readRDS)


# Get top 500 bins for every topic  ---------------------------------------

tm.result.lst <- lapply(lda.objs, function(x) posterior(x))

# Get top bins  -----------------------------------------------------------

top.features.common.lst <- lapply(jmarks.single, function(jmarktmp){
  tm.result <- tm.result.lst[[jmarktmp]]
  cnames <- colnames(tm.result$terms)
  top.features.lst <- apply(tm.result$terms, 1, function(jrow) GetTopFeatureNames(jrow, cnames, topbins.keep))
  top.features.common <- unique(unlist(top.features.lst))
  return(top.features.common)
})


features.keep.merged <- unlist(top.features.common.lst)

# Write new count mat, filtered by features  ------------------------------

count.mat.filt.lst <- lapply(jmarks, function(jmark){
  count.mat.tmp <- count.mats[[jmark]][features.keep.merged, ]
})


# Write output  -----------------------------------------------------------

# count.mats.filt and bedfiles of all the locations

for (jmarktmp in jmarks){
  outmat.tmp <- file.path(outdir, paste0("countmat_featurefilt.", jmarktmp, ".", Sys.Date(), ".rds"))
  # write mats
  saveRDS(count.mat.filt.lst[[jmarktmp]], file = outmat.tmp)
}
