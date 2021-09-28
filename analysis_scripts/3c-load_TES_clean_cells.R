# Jake Yeung
# Date of Creation: 2021-07-09
# File: ~/projects/scChIX/analysis_scripts/3c-load_TES_clean_cells.R
#


rm(lis=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)




# Load TES  ---------------------------------------------------------------

# jmarks <- c("K36", "K9m3")
# names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_TES/lda_outputs.TES_counts.K36.2021-06-30.K-30.binarize.FALSE/ldaOut.TES_counts.K36.2021-06-30.K-30.Robj")
load(inf, v=T)


# Load metas  -------------------------------------------------------------


inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_NN_15_check_plates/dbl_k36_k9m3_cleaned/celltyping_output_filt.K36.2021-07-07.rds")
dat.meta <- readRDS(inf.meta)
cells.keep <- dat.meta$cell

cols.keep <- colnames(count.mat) %in% cells.keep
count.mat.filt <- count.mat[, cols.keep]

dim(count.mat)
dim(count.mat.filt)


# SAave

outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_NN_15_check_plates/TES_k36_cleaned")
outrds <- file.path(outdir, paste0("countmat_TES_cleaned.K36.2021-07-09.rds"))
saveRDS(count.mat.filt, outrds)
