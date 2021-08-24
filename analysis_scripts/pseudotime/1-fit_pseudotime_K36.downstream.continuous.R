# Jake Yeung
# Date of Creation: 2021-08-07
# File: ~/projects/scChIX/analysis_scripts/pseudotime/6-fit_pseudotime_K36.downstream.continuous.R
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

library(scchicFuncs)
library(scChIX)
library(JFuncs)

# Functions ---------------------------------------------------------------




# Constants ---------------------------------------------------------------

jmark <- "K36"
outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs"
jprefix <- "var_filtered"
jsuffix <- "manual2nocenter_K36_K9m3_K36-K9m3"
jname <- paste(jprefix, jsuffix, sep = "_")
outdir <- file.path(outmain, jsuffix)

jdate <- "2021-08-07"

# outf <- file.path(outdir, paste0("glm_poisson_fits_output.permute.", jsuffix, ".", jmark, ".", jdate, ".RData"))
outf <- file.path(outdir, paste0("glm_poisson_fits_output.permute.seed_1.", jsuffix, ".", jmark, ".", jdate, ".RData"))
assertthat::assert_that(file.exists(outf))


# check output
load(outf, v=T)

jrowname <- names(jfits.lst)
names(jrowname) <- jrowname

pval.lst <- lapply(jfits.lst, function(jfit){
  jfit$pval
})

slope.lst <- lapply(jfits.lst, function(jfit){
  jfit$ptime.Estimate
})

dat.slope.pval <- lapply(jrowname, function(jrowname){
  dat.out <- data.frame(bin = jrowname, pval = pval.lst[[jrowname]], slope = slope.lst[[jrowname]], stringsAsFactors = FALSE)
}) %>%
  bind_rows()

ggplot(dat.slope.pval, aes(x = slope)) + geom_histogram(bins = 1000) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 2, color = "red")

ggplot(dat.slope.pval, aes(x = pval)) + geom_histogram(bins = 100) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(0, 1))

ggplot(dat.slope.pval, aes(x = slope, y = -log10(pval))) + geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 1, color = "red") +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1, color = "red") +
  coord_cartesian(xlim = c(-3, 3))

dat.slope.pval.sorted <- dat.slope.pval %>%
  arrange(pval) %>%
  mutate(qval = p.adjust(pval)) %>%
  filter(slope < 0)


print(head(dat.slope.pval.sorted))

jbin <- "chr3:106150000-106200000"
jbin <- "chr9:3000000-3050000"
jbin <- dat.slope.pval.sorted$bin[5]

jbin <- dat.slope.pval.sorted$bin[3]
jbin <- dat.slope.pval.sorted$bin[4]
jbin <- dat.slope.pval.sorted$bin[5]
jbin <- dat.slope.pval.sorted$bin[6]
jbin <- dat.slope.pval.sorted$bin[7]

jbin <- dat.slope.pval.sorted$bin[2]
jbin <- dat.slope.pval.sorted$bin[3]

jbin <- dat.slope.pval.sorted$bin[2]
jbin <- dat.slope.pval.sorted$bin[3]
jbin <- dat.slope.pval.sorted$bin[4]


jbin <- dat.slope.pval.sorted$bin[1]

jcheck <- data.frame(counts = count.mat.filt[jbin, ], cell = colnames(count.mat.filt), stringsAsFactors = FALSE) %>%
  left_join(., dat.annots.filt)  %>%
  left_join(., subset(dat.umap, select = c(cell, umap1, umap2))) %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]]) %>%
  ungroup() %>%
  mutate(ptime.factor = factor(ptime, levels = c(9.5, 10.5, 11.5)))

jcheck.sum <- jcheck %>%
  group_by(ptime.factor) %>%
  summarise(nbr.nonzeros = nnzero(counts),
            ncells = length(counts)) %>%
  ungroup() %>%
  mutate(frac.nonzeros = nbr.nonzeros / ncells)

ggplot(jcheck, aes(x = ptime.factor, y = log(counts + 1))) +
  geom_point(alpha = 0.1) +
  geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jcheck.sum, aes(x = ptime.factor, y = frac.nonzeros)) +
  geom_col() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot on UMAP

ggplot(jcheck, aes(x = umap1, y = umap2, color = counts)) +
  theme_bw() +
  facet_wrap(~ptime.factor) +
  geom_point() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
