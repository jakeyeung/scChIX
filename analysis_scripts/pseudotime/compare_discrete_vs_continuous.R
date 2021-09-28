# Jake Yeung
# Date of Creation: 2021-08-05
# File: ~/projects/scChIX/analysis_scripts/pseudotime/compare_discrete_vs_continuous.R
#

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load objs ---------------------------------------------------------------

inf1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/manual2nocenter_K36_K9m3_K36-K9m3/glm_poisson_fits_output.manual2nocenter_K36_K9m3_K36-K9m3.2021-08-05.RData"
load(inf1, v=T)

jfits.lst1 <- jfits.lst
dat.annots.filt1 <- dat.annots.filt
ncuts.cells1 <- ncuts.cells
count.mat.filt1 <- count.mat.filt
dat.umap.filt1 <- dat.umap.filter

inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/manual2nocenter_K36_K9m3_K36-K9m3/glm_poisson_fits_output.discrete.manual2nocenter_K36_K9m3_K36-K9m3.2021-08-05.RData"
load(inf2, v=T)

jfits.lst2 <- jfits.lst
dat.annots.filt2 <- dat.annots.filt
ncuts.cells2 <- ncuts.cells
count.mat.filt2 <- count.mat.filt
dat.umap.filt2 <- dat.umap.filter


# are significant hits in continuous same as in discrete?
jnames <- names(jfits.lst)
names(jnames) <- jnames

pvals.vec1 <- unlist(lapply(jfits.lst1, function(x) x$pval))
pvals.vec2 <- unlist(lapply(jfits.lst2, function(x) x$pval))

pvals.vec1 <- unlist(lapply(jfits.lst1, function(x) x$ptime.Estimate))
pvals.vec2 <- unlist(lapply(jfits.lst2, function(x) x$ptime.Estimate / 2))

plot(pvals.vec1, pvals.vec2, xlim = c(-2, 2), ylim = c(-2, 2), pch = 20)


# jcheck <- data.frame(counts = count.mat.filt[jbin, ], cell = colnames(count.mat.filt), stringsAsFactors = FALSE) %>%
#   left_join(., dat.annots.filt)  %>%
#   left_join(., subset(dat.umap, select = c(cell, umap1, umap2))) %>%
#   rowwise() %>%
#   mutate(stage = strsplit(cell, split = "-")[[1]][[1]]) %>%
#   ungroup() %>%
#   mutate(ptime.factor = factor(ptime, levels = c(9.5, 10.5, 11.5)))
#
# jcheck.sum <- jcheck %>%
#   group_by(ptime.factor) %>%
#   summarise(nbr.nonzeros = nnzero(counts),
#             ncells = length(counts)) %>%
#   ungroup() %>%
#   mutate(frac.nonzeros = nbr.nonzeros / ncells)


# Write metafile for plotting hits  ---------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_pseudotime"
outf <- file.path(outdir, paste0("metadata_pseudotime.", Sys.Date(), ".txt"))

dat.umap.ordered.final <- dat.umap.filt1 %>%
  left_join(., dat.annots.filt1) %>%
  arrange(ptime) %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]]) %>%
  mutate(stage.numeric = as.numeric(gsub("p", ".", gsub("^E", "", stage))))


fwrite(dat.umap.ordered.final, file = outf, sep = "\t")




