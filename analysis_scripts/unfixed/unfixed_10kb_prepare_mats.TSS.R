# Jake Yeung
# Date of Creation: 2021-09-06
# File: ~/projects/scChIX/analysis_scripts/unfixed/unfixed_10kb_prepare_mats.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("K4m1", "K27m3", "K4m1_K27m3"); names(jmarks) <- jmarks

# Load 10kb mats ---------------------------------------------------------

infs <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/2021-03-17_redo_count_tables.TSS_allgenes/countTables.unfixed.TSS.refseq/all_BM_", jmark, "_200119.mq_40.TSS.winsize_10000.csv.gz"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})


dats <- lapply(infs, function(jinf){
  dat <- ReadMatTSSFormat(jinf)
})


# Filter cells ------------------------------------------------------------

inf.rdata <- file.path(hubprefix, "jyeung/data/dblchic/from_rstudio/primetime/unfixed_louvain2/BM_UnfixedLouvain2.FinalCellClusterTable.2020-03-21.RData")
load(inf.rdata, v=T)

cells.keep <- unique(dat.final.annots$cell)

dats.cellfilt <- lapply(jmarks, function(jmark){
  print(jmark)
  mat.tmp <- dats[[jmark]]
  cols.keep <- colnames(dats[[jmark]]) %in% cells.keep
  mat.tmp.filt <- mat.tmp[, cols.keep]
  return(mat.tmp.filt)
})

lapply(dats, dim)
lapply(dats.cellfilt, dim)


# Filter bins (optional?) -------------------------------------------------


# jcutoff.log <- 6
jcutoff.log <- 5.5
top.genes <- 5000

rows.keep.lst <- lapply(jmarks, function(jmark.tmp){
  print(jmark.tmp)
  # jmark.tmp <- jmarks[[3]]
  mat.tmp <- dats.cellfilt[[jmark.tmp]]

  nvec <- colSums(mat.tmp)

  gdevs <- apply(mat.tmp, 1, function(xvec){
    scchicFuncs::binomial_deviance(x = xvec, p = sum(xvec) / sum(nvec), n = nvec)
  })

  gdevs.sorted <- sort(gdevs, decreasing = TRUE)

  rows.keep.tmp <- names(gdevs.sorted[1:top.genes])

  # plot(density(log(gdevs)), main = jmark.tmp)
  # abline(v = jcutoff.log, col = 'blue', lty = 'dotted')
  # (jcutoff <- exp(jcutoff.log))
#
  # rows.keep.tmp <- gdevs[which(gdevs > jcutoff)]
  return(rows.keep.tmp)
})

print(lapply(rows.keep.lst, length))

# common rows
rows.keep.common <- unique(unlist(rows.keep.lst))

print(length(rows.keep.common))


# Output ------------------------------------------------------------------

dats.cellfilt.binfilt <- lapply(dats.cellfilt, function(jdat){
  jdat[rows.keep.common, ]
  cells.keep <- colSums(jdat.filt) > 0
  return(jdat.filt[, cells.keep])
})

lapply(dats.cellfilt.binfilt, dim)
lapply(dats.cellfilt, dim)

# write outputs
outdir.gw <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/count_mats_for_LDA.10kb/10kb_TSS")
outdir.gwfilt <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/count_mats_for_LDA.10kb/10kb_TSS_filt")
dir.create(outdir.gw)
dir.create(outdir.gwfilt)


for (jmark in jmarks){
  print(jmark)
  outf.gw <- file.path(outdir.gw, paste0("count_mat.10kb_TSS.", jmark, ".rds"))
  outf.gwfilt <- file.path(outdir.gwfilt, paste0("count_mat.10kb_TSS_filt.", jmark, ".rds"))
  saveRDS(dats.cellfilt[[jmark]], file = outf.gw)
  saveRDS(dats.cellfilt.binfilt[[jmark]], file = outf.gwfilt)
}

