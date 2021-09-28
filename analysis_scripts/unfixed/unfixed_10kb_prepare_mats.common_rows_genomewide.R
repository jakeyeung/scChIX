# Jake Yeung
# Date of Creation: 2021-09-12
# File: ~/projects/scChIX/analysis_scripts/unfixed/unfixed_10kb_prepare_mats.common_rows_genomewide.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("K4m1", "K27m3", "K4m1-K27m3"); names(jmarks) <- jmarks

# Load mats ---------------------------------------------------------------


indir <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/count_mats_for_LDA.10kb/10kb_genomewide")

inf.mats <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, paste0("count_mat.10kb_genomewide.", jmark, ".rds"))
  return(inf.tmp)
})

mats <- lapply(inf.mats, function(jinf){
  readRDS(jinf)
})


rows.lst <- lapply(mats, function(x) rownames(x))
common.rows <- Reduce(intersect, rows.lst)

mats.commonrows <- lapply(mats, function(jmat){
  jmat[common.rows, ]
})

rnames.check.lst <- lapply(mats.commonrows, rownames)
assertthat::assert_that(length(Reduce(intersect, rnames.check.lst)) == length(common.rows))

# rnames.lst <- list(colnames(tm.result.lst[[jmark1]]$terms), colnames(tm.result.lst[[jmark2]]$terms), rownames(check3))
# rnames.common <- Reduce(f = intersect, x = rnames.lst)
#
# jbin <- "chr12:19530000-19540000"
# mats.commonrows[[1]][jbin, ]
# mats.commonrows[[2]][jbin, ]
# mats.commonrows[[3]][jbin, ]

# Check -------------------------------------------------------------------




# Save outputs ------------------------------------------------------------

for(jmark in jmarks){
  out.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_inputs/countmats/countmat_var_filt.", jmark, ".rds"))
  print(dim(mats.commonrows[[jmark]]))
  saveRDS(object = mats.commonrows[[jmark]], file = out.tmp)
}
