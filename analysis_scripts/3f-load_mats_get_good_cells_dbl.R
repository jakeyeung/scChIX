# Jake Yeung
# Date of Creation: 2021-07-11
# File: ~/projects/scChIX/analysis_scripts/3f-load_mats_get_good_cells_dbl.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks

# Load metas  -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
infs.meta <- lapply(jmarks, function(jmark){
  inf.meta.tmp <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_LDA/dbl_cleaned/celltyping_output_filt.", jmark, ".2021-07-10.rds"))
  assertthat::assert_that(file.exists(inf.meta.tmp))
  return(inf.meta.tmp)
})

dats.meta <- lapply(infs.meta, function(jinf){
  readRDS(jinf)
})

print(lapply(dats.meta, dim))


# Load data  --------------------------------------------------------------

# get K36-K9m3

inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/tagged_bams/counts_tables_topfeatures_K36_genebodies_K9m3_bins")
fnames.all <- list.files(inmain, pattern = "*.txt")

fnames.dbl.i <- sapply(fnames.all, function(f) grepl(pattern = "K36-K9", f))
fnames.dbl <- fnames.all[fnames.dbl.i]

mats.lst.dbl <- sapply(fnames.dbl, function(fname){
  inf.tmp <- file.path(inmain, fname)
  ReadMatTSSFormat(inf.tmp)
})

# get K36s
fnames.k36.i <- sapply(fnames.all, function(f) grepl(pattern = "^K36", strsplit(f, split = "-")[[1]][[3]]) & grepl(pattern = "^19", strsplit(f, split = "-")[[1]][[4]]))
fnames.k36 <- fnames.all[fnames.k36.i]

fnames.k9m3.i <- sapply(fnames.all, function(f) grepl(pattern = "^K9", strsplit(f, split = "-")[[1]][[3]]) & grepl(pattern = "^19", strsplit(f, split = "-")[[1]][[4]]))
fnames.k9m3 <- fnames.all[fnames.k9m3.i]

mats.lst.k36 <- sapply(fnames.k36, function(fname){
  inf.tmp <- file.path(inmain, fname)
  ReadMatTSSFormat(inf.tmp)
})

mats.lst.k9m3 <- sapply(fnames.k9m3, function(fname){
  inf.tmp <- file.path(inmain, fname)
  ReadMatTSSFormat(inf.tmp)
})

mats.lst.both <- c(mats.lst.k36, mats.lst.k9m3, mats.lst.dbl)

all.rnames.both.lst <- lapply(mats.lst.both, function(jmat) rownames(jmat))
all.rnames.both.merged <- sort(unique(unlist(all.rnames.both.lst)))

mat.merged.both <- cbind.fill.lst(mats.lst.both, all.rnames = all.rnames.both.merged, fill = 0)



# Get cells split into two mats -------------------------------------------

cells.keep.lst <- lapply(dats.meta, function(jdat){
  jdat$cell
})

mats.merged.lst <- lapply(jmarks, function(jmark){
  cells.keep.tmp <- cells.keep.lst[[jmark]]
  mat.merged.both[, cells.keep.tmp]
})

print(lapply(mats.merged.lst, dim))

print(lapply(cells.keep.lst, length))


# Write mats  -------------------------------------------------------------

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/coords_filtered/cells_filtered"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/coords_filtered_K36_K9m3_merged_with_dbl"
dir.create(outdir)

outrds.lst <- lapply(jmarks, function(jmark){
  outrds.tmp <- file.path(outdir, paste0("countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.", Sys.Date(), ".", jmark, ".rds"))
  saveRDS(mats.merged.lst[[jmark]], file = outrds.tmp)
})
# outrds1 <- file.path(outdir, paste0("countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.", Sys.Date(), ".", jmarks[[1]], ".rds"))
# outrds2 <- file.path(outdir, paste0("countmat_featuresfilt_cellsfilt_K36_genebodies_K9m3_bins.", Sys.Date(), ".", jmarks[[2]], ".rds"))
# saveRDS(mats.merged.k9m3, file = outrds2)
