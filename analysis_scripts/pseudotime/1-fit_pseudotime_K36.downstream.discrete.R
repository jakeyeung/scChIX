# Jake Yeung
# Date of Creation: 2021-08-07
# File: ~/projects/scChIX/analysis_scripts/pseudotime/6-fit_pseudotime_K36.downstream.discrete.R
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

# outf <- file.path(outdir, paste0("glm_poisson_fits_output.discrete.", jsuffix, ".", jmark, ".", jdate, ".RData"))
# assertthat::assert_that(file.exists(outf))

outf <- file.path(outdir, paste0("glm_poisson_fits_output.discrete.permute.", jsuffix, ".", jmark, ".", jdate, ".RData"))
assertthat::assert_that(file.exists(outf))

# check output
load(outf, v=T)

print(head(dat.umap.filt))


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



# Plot raw  ---------------------------------------------------------------



