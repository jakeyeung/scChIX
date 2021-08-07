# Jake Yeung
# Date of Creation: 2021-08-05
# File: ~/projects/scChIX/analysis_scripts/pseudotime/filter_mats_for_LDA.R
# Filter out islands in K9me3 and then make mat for LD

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load mat ----------------------------------------------------------------

out.check <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/manual2nocenter_K36_K9m3_K36-K9m3/glm_poisson_fits_output.manual2nocenter_K36_K9m3_K36-K9m3.2021-08-04.RData"
load(out.check, v=T)


# Save rds  ---------------------------------------------------------------

outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_count_tables_filtered/countmat_K9me3_pseudotime_filt.", Sys.Date(), ".rds")
saveRDS(count.mat.filt, file = outf)
