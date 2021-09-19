# Jake Yeung
# Date of Creation: 2021-08-11
# File: ~/projects/scChIX/analysis_scripts/3f-load_10kb_filter_cells.R
# Load 10kb bin raw counts, filter cells, write count tables for LDA

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)



# Load 10kb  --------------------------------------------------------------

jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/count_tables_dedup_from_merged/counts_tables_10000")

# load tables
infs.lst <- lapply(jmarks, function(jmark){
  fname.mark <- paste0("gastru_merged_", jmark, ".rowsdeduped.tagged.sorted.countTable.binsize_10000.csv")
  inf.mark <- file.path(indir, fname.mark)
  assertthat::assert_that(file.exists(inf.mark))
  return(inf.mark)
})

mats.lst <- lapply(infs.lst, function(jinf){
  print(jinf)
  ReadMatSlideWinFormat(jinf)
})

# get common rows
rows.lst <- lapply(mats.lst, function(jmat){
  rownames(jmat)
})

rows.common <- Reduce(intersect, x = rows.lst)


# Filter cells ------------------------------------------------------------

# load cells

jname <- "var_filtered_manual2nocenternoE8unifyK36_K36_K27_K36-K27"
# inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenternoE8unifyK36_K36_K27_K36-K27/celltyping_K36_first_try.2021-08-02.txt")
inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/", jname, "/metamerged_bothmarks_K36_first_try.2021-08-11.txt")
dat.meta <- fread(inf.meta)

cells.keep <- dat.meta$cell

mats.rowscommon.cellsfilt.lst <- lapply(mats.lst, function(jmat){
  cells.filt <- colnames(jmat)%in% cells.keep
  return(jmat[rows.common, cells.filt])
})

# Write mat outputs ------------------------------------------------------

# jname <- "K36_K27_K36-K27"
outmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/binsize_10000")
dir.create(outmain)
outdir <- file.path(outmain, jname)
dir.create(outdir)

metatmp <- paste0("metadata_var_filt.K36_K27_K36-K27.", Sys.Date(), ".rds")
metaftmp <- file.path(outdir, metatmp)
saveRDS(object = dat.meta, file = metaftmp)

for (jmark in jmarks){
  print(jmark)
  fnametmp <- paste0("countmat_var_filt.", jmark, ".", Sys.Date(), ".rds")
  outftmp <- file.path(outdir, fnametmp)
  mattmp <- mats.rowscommon.cellsfilt.lst[[jmark]]

  print(dim(mattmp))

  saveRDS(object = mattmp, file = outftmp)
}

