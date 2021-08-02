# Jake Yeung
# Date of Creation: 2021-06-30
# File: ~/projects/scChIX/analysis_scripts/1-load_bins10kb_mats_filter_good_cells.R
#


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

# jsuffix <- "TSS10kb"
jsuffix <- "10000"
hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmarks <- c("K4m1", "K36", "K27"); names(jmarks) <- jmarks
jmarks <- c("K36", "K27", "K9m3", "K36-K9m3", "K36-K27"); names(jmarks) <- jmarks

inmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables/counts_tables_", jsuffix))
assertthat::assert_that(dir.exists(inmain))

# Load mat  ---------------------------------------------------------------

infs.rds <- lapply(jmarks, function(jmark){
  # inf.tmp <- file.path(inmain, paste0(jmark, "/cbind_out/", jmark, ".counts_tables_", jsuffix, ".rds"))
  inf.tmp <- file.path(inmain, paste0(jmark, "/cbind_out/", jmark, ".countTable.binsize_", jsuffix, ".rds"))
  print(inf.tmp)
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

mats.lst <- lapply(infs.rds, function(inf){
  readRDS(inf)
})


# Filt cells --------------------------------------------------------------

# inf.meta.k36 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables/K36_K27_K36-K27/K36_K27_K36-K27.meta_data.50000.K27.2021-06-28.txt")

indir.meta.k36 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_filtered_counts_varcutoffmin_0.03/K36_K9m3_K36-K9m3")
indir.meta.k27 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_filtered_counts_varcutoffmin_0.03/K36_K27_K36-K27")

infs.meta.k36 <- list.files(indir.meta.k36, pattern = "*.txt", full.names = TRUE)
infs.meta.k27.all <- list.files(indir.meta.k27, pattern = "*.txt", full.names = TRUE)
infs.meta.k27.keep <- sapply(infs.meta.k27, function(x) strsplit(x = basename(x), split = "\\.")[[1]][[4]] != "K36")
infs.meta.k27 <- infs.meta.k27.all[infs.meta.k27.keep]
infs.meta.all <- c(infs.meta.k36, infs.meta.k27)


dats.meta <- lapply(infs.meta.all, fread) %>%
  bind_rows()

good.cells <- dats.meta$cell


mats.filt.lst <- lapply(mats.lst, function(jmat){
  cells.keep <- colnames(jmat) %in% good.cells
  jmat.filt <- jmat[, cells.keep]
  print(dim(jmat))
  print(dim(jmat.filt))
  return(jmat.filt)
})


# Output  -----------------------------------------------------------------


outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/binsize_", jsuffix, "_filtered_count_tables_filtered_counts_varcutoffmin_0.03"))
dir.create(outdir)
for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/binsize_", jsuffix, "_filtered_count_tables_filtered_counts_varcutoffmin_0.03/binsize_", jsuffix, ".", jmark, ".", Sys.Date(), ".rds"))
  jmat.tmp <- mats.filt.lst[[jmark]]
  print(dim(jmat.tmp))
  saveRDS(jmat.tmp, file = outrds)
}

