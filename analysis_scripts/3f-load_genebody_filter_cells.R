# Jake Yeung
# Date of Creation: 2021-08-13
# File: ~/projects/scChIX/analysis_scripts/3f-load_genebody_filter_cells.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)
library(scchicFuncs)

jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"
inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/tagged_bams")
# jbinsize <- "genebody50kbmax"
# jbinsize <- "20000"
# jbinsize <- "30000"
# jbinsize <- "40000"
jbinsize <- "10000"
indir.counts <- file.path(inmain, paste0("counts_tables_", jbinsize))
dir.exists(indir.counts)

jprefix <- "var_filtered_manual2nocenternoE8unifyK36"
jstr <- "K36_K27_K36-K27"
jname <- paste(jprefix, jstr, sep = "_")

# Load plate and good cells  ----------------------------------------------

indir.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/", jprefix, "/K36_K27_K36-K27"))
dir.exists(indir.meta)

infs.meta <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir.meta, paste0("meta_var_filt.", jmark, ".2021-07-29.rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dats.meta <- lapply(infs.meta, function(inf){
  readRDS(inf)
})

cells.keep <- unique(unlist(lapply(dats.meta, function(jdat) jdat$cell)))


# Get experis per mark  ---------------------------------------------------

experi.lst <- lapply(dats.meta, function(jdat){
  unique(sapply(jdat$cell, function(x) ClipLast(x = x, jsep = "_")))
})

# Load mats, filt cells ---------------------------------------------------


# indir <- file.path(inmain, "")

mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  experi.vec <- experi.lst[[jmark]]
  names(experi.vec) <- experi.vec
  mats.byplate.lst <- lapply(experi.vec, function(experi){
    print(experi)
    if (jbinsize == "genebody50kbmax"){
      fpath <- file.path(indir.counts, paste0(experi, ".sorted.tagged.bam.count_table_", jbinsize, ".txt"))
      mat.tmp <- scchicFuncs::ReadMatTSSFormat(fpath)
    } else {
      fpath <- file.path(indir.counts, paste0(experi, ".sorted.tagged.countTable.binsize_", jbinsize, ".csv"))
      mat.tmp <- scchicFuncs::ReadMatSlideWinFormat(fpath)
    }
  })
  # bind by marks
  all.rnames <- lapply(mats.byplate.lst, function(jmat) rownames(jmat))
  all.rnames <- unique(unlist(all.rnames))
  mats.byplate.merged <- cbind.fill.lst(mats.lst = mats.byplate.lst, all.rnames = all.rnames, fill = 0)
})


# get common rows
rows.lst <- lapply(mats.lst, function(jmat){
  rownames(jmat)
})

rows.common <- Reduce(intersect, x = rows.lst)


mats.rowscommon.cellsfilt.lst <- lapply(mats.lst, function(jmat){
  cells.filt <- colnames(jmat)%in% cells.keep
  return(jmat[rows.common, cells.filt])
})

# Write mat outputs ------------------------------------------------------

outmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/binsize_", jbinsize))
dir.create(outmain)
outdir <- file.path(outmain, jname)
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  fnametmp <- paste0("countmat_var_filt.", jmark, ".", Sys.Date(), ".rds")
  outftmp <- file.path(outdir, fnametmp)
  mattmp <- mats.rowscommon.cellsfilt.lst[[jmark]]
  print(dim(mattmp))
  saveRDS(object = mattmp, file = outftmp)
}

