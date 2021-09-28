# Jake Yeung
# Date of Creation: 2021-08-09
# File: ~/projects/scChIX/analysis_scripts/6-make_metas_for_plotting.R
# Make metas with proper column names for plotting

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)




# Load metas --------------------------------------------------------------

jmarks <- c("K9m3", "K36"); names(jmarks) <- jmarks


hubprefix <- "/home/jyeung/hub_oudenaarden"

infs <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/celltyping_", jmark, "_first_try.2021-08-02.txt"))
  assertthat::assert_that(file.exists(inf.meta))
  return(inf.meta)
})

dat.metas <- lapply(infs, function(inf) fread(inf))

dat.metas[[1]]$celltype <- dat.metas[[1]]$cluster

# Add colors --------------------------------------------------------------

# "clustercol"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

dat.metas <- lapply(dat.metas, function(jdat){
  jdat$celltype.factor <- factor(jdat$celltype)
  # jdat$cluster <- as.character(jdat$celltype)
  jdat$celltype <- NULL
  jdat$cluster <- factor(as.character(jdat$stage), levels = c("E9p5", "E10", "E10p5", "E11p5"))
  jdat$clustercol <- sapply(as.numeric(jdat$cluster), function(i) cbPalette[[i]])
  return(jdat)
})

ggplot(dat.metas$K9m3, aes(x = umap1, y = umap2, color = clustercol)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.metas$K36, aes(x = umap1, y = umap2, color = clustercol)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Load raw counts to calculate total cuts  --------------------------------

inmain.raw <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3")
infs.raw <- lapply(jmarks, function(jmark){
  jsuffix <- paste0("lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj")
  inf.raw <- file.path(inmain.raw, jsuffix)
  assertthat::assert_that(file.exists(inf.raw))
  return(inf.raw)
})

count.raw.lst <- lapply(infs.raw, function(jinf){
  load(jinf, v=T)
  return(count.mat)
})

cuts.total.lst <- lapply(count.raw.lst, function(jcount){
  data.frame(cell = colnames(jcount), cuts_total = colSums(jcount), stringsAsFactors = FALSE)
})

dat.metas <- lapply(jmarks, function(jmark){
  left_join(dat.metas[[jmark]], cuts.total.lst[[jmark]])
})

# sort by pseudotime for K9
inf.pseudo <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_pseudotime/metadata_pseudotime.2021-08-05.txt")
dat.pseudo <- fread(inf.pseudo)

dat.metas$K9m3 <- left_join(dat.metas$K9m3, subset(dat.pseudo, select = c(cell, ptime))) %>%
  arrange(ptime) %>%
  filter(!is.na(ptime))

dat.metas$K36 <- dat.metas$K36 %>%
  arrange(cluster, type)


# Write outputs -----------------------------------------------------------

outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots")
for (jmark in jmarks){
  outpdf <- file.path(outdir, paste0("clusterplot_pseudotime_K36-K9me3.", jmark, ".", Sys.Date(), ".pdf"))
  outtxt <- file.path(outdir, paste0("metadata_pseudotime_K36-K9me3.", jmark, ".", Sys.Date(), ".txt"))

  pdf(outpdf, useDingbats = FALSE)
  m1 <- ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = clustercol)) +
    geom_point() +
    ggtitle(jmark) +
    theme_bw() +
    scale_color_identity() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)

  fwrite(dat.metas[[jmark]], file = outtxt, sep = "\t")

  dev.off()
}


