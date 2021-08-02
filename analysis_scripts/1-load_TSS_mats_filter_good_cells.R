# Jake Yeung
# Date of Creation: 2021-06-30
# File: ~/projects/scChIX/analysis_scripts/1-load_TSS_mats_filter_good_cells.R
#


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jsuffix <- "TSS10kb"
hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("K4m1", "K36", "K27"); names(jmarks) <- jmarks
jmarks <- c("K36", "K27"); names(jmarks) <- jmarks

inmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables/counts_tables_", jsuffix))
assertthat::assert_that(dir.exists(inmain))

# Load mat  ---------------------------------------------------------------

infs.rds <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(inmain, paste0(jmark, "/cbind_out/", jmark, ".counts_tables_", jsuffix, ".rds"))
  print(inf.tmp)
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

mats.lst <- lapply(infs.rds, function(inf){
  readRDS(inf)
})


# Filt cells --------------------------------------------------------------

# inf.meta.k36 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables/K36_K27_K36-K27/K36_K27_K36-K27.meta_data.50000.K27.2021-06-28.txt")
inf.meta.k36 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_filtered_counts_varcutoffmin_0.03/K36_K9m3_K36-K9m3/meta_data.50000.K36.2021-06-29.txt")
inf.meta.k27 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_filtered_counts_varcutoffmin_0.03/K36_K27_K36-K27/K36_K27_K36-K27.meta_data.50000.K27.2021-06-30.txt")

dat.meta.k36 <- fread(inf.meta.k36)
dat.meta.k27 <- fread(inf.meta.k27)

cells.k36 <- dat.meta.k36$cell
cells.k27 <- dat.meta.k27$cell

good.cells <- c(cells.k36, cells.k27)

mats.filt.lst <- lapply(mats.lst, function(jmat){
  cells.keep <- colnames(jmat) %in% good.cells
  jmat.filt <- jmat[, cells.keep]
  print(dim(jmat))
  print(dim(jmat.filt))
  return(jmat.filt)
})


# Output  -----------------------------------------------------------------


outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/", jsuffix, "_filtered_count_tables_filtered_counts_varcutoffmin_0.03"))
dir.create(outdir)
for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/", jsuffix, "_filtered_count_tables_filtered_counts_varcutoffmin_0.03/", jsuffix, ".", jmark, ".", Sys.Date(), ".rds"))
  jmat.tmp <- mats.filt.lst[[jmark]]
  print(dim(jmat.tmp))
  saveRDS(jmat.tmp, file = outrds)
}

