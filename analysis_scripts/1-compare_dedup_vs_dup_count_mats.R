# Jake Yeung
# Date of Creation: 2021-07-31
# File: ~/projects/scChIX/analysis_scripts/1-compare_dedup_vs_dup_count_mats.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

# Compare dup vs dedup'd count mats ---------------------------------------

jprefix <- "/home/jyeung/hub_oudenaarden"

# dup countmat is already in nice .rds format
inf.dup <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/countmat_output_filt.K36-K9m3.2021-07-23.rds"
mat.dup <- readRDS(inf.dup)

inf.dedup <- file.path(jprefix, "jyeung/data/dblchic/gastrulation/count_tables_dedup_from_merged/counts_tables_50000/gastru_merged_K36-K9m3.rowsdeduped.tagged.sorted.countTable.binsize_50000.csv")
mat.dedup <- ReadMatSlideWinFormat(inf.dedup)

cells.keep <- colnames(mat.dup)
cols.keep <- colnames(mat.dedup) %in% cells.keep

mat.dedup.filt <- mat.dedup[, cols.keep]

i <- sample(cells.keep, size = 1)
print(i)
table(mat.dedup.filt[, i])
table(mat.dup[, i])


plot(density(log10(as.matrix(mat.dedup.filt[, i]) + 1)), main = i)
plot(density(log10(mat.dup[, i] + 1)), main = i)

