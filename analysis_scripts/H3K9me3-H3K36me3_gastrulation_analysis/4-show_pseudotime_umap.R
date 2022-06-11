# Jake Yeung
# Date of Creation: 2021-10-04
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/4-show_pseudotime_umap.R
# 



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load meta ---------------------------------------------------------------

hubprefix <- "/Users/yeung/hub_oudenaarden"
inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/metadata_cleaned/meatdata_cleaned.2021-10-04.rds"))
assertthat::assert_that(file.exists(inf.meta))

dat.meta <- readRDS(inf.meta)


# Load raw counts H3K9me3  ------------------------------------------------

inf.mat <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/scchix_unmixing_downstream/scchix_inputs_clstr_by_celltype-merged_mat.K9m3.rds")
mat <- readRDS(inf.mat)


# Show timepoints for H3K9me3  --------------------------------------------

# filter out
clsts.remove <- c("Erythroid", "WhiteBloodCells")

dat.filt <- subset(dat.meta, mark == "K9m3" & !cluster %in% clsts.remove)

ggplot(dat.filt, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cells.keep <- dat.filt$cell

mat.filt <- mat[, cells.keep]


# Write output ------------------------------------------------------------


outf <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/downstream/count_mats_rds_objs/count_mat_gastru_H3K9me3_no_blood.2021-10-04.rds")
saveRDS(mat.filt, file = outf)


