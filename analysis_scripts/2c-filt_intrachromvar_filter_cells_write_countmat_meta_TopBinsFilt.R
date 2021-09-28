# Jake Yeung
# Date of Creation: 2021-07-15
# File: ~/projects/scChIX/analysis_scripts/2c-filt_intrachromvar_filter_cells_write_countmat_meta_TopBinsFilt.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# jquantile <- 0.15
jquantile <- 0.3


# Load  -------------------------------------------------------------------

jdate.output <- "2021-07-15"
jdate <- "2021-07-14"
# jmarks <- c("K36", "K27", "K36-K27")
jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

hubprefix <- "/home/jyeung/hub_oudenaarden"
inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/", jmarks[[2]], "/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt")
dir.exists(inmain)

inmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/var_filtered_", jquantile), jstr)
inmain.counts <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/coords_filtered"), jmarks[[2]])
outmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/coords_var_filtered_", jquantile))
dir.create(outmain)
outdir <- file.path(outmain, jstr)
dir.create(outdir)


# Load metas  -------------------------------------------------------------

infs.meta <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(inmain, paste0("celltyping_var_filt.", jmark, ".2021-07-15.rds"))
  assertthat::assert_that(file.exists(inf.meta))
  return(inf.meta)
})


# Load counts from TopBinsFilt --------------------------------------------


infs.countmat <- lapply(jmarks, function(jmark){
  inf.countmat <- file.path(inmain.counts, paste0("countmat_featurefilt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.countmat))
  return(inf.countmat)
})

dat.umap.lst <- lapply(infs.meta, readRDS)
countmats.lst <- lapply(infs.countmat, readRDS)


for (jmark in jmarks){
  print(jmark)
  # write dat impute vars

  # write filtered count mat
  cells.keep.tmp <- dat.umap.lst[[jmark]]$cell
  dat.countmat.tmp <- countmats.lst[[jmark]][, cells.keep.tmp]
  outf.countmat.tmp <- file.path(outdir, paste0("countmat_var_filt_TopBinsFilt.", jmark, ".", jdate.output, ".rds"))

  print(dim(countmats.lst[[jmark]]))
  print(dim(dat.countmat.tmp))

  assertthat::assert_that(ncol(dat.countmat.tmp) < ncol(countmats.lst[[jmark]]))

  saveRDS(dat.countmat.tmp, outf.countmat.tmp)
}

