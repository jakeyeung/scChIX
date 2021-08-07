# Jake Yeung
# Date of Creation: 2021-08-04
# File: ~/projects/scChIX/analysis_scripts/6-fit_pseudotime_K9me3.R
# Fit pseudotime on single+dbl UMAP


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


FitGlmRowPtime.withse <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE, with.se = FALSE){
  # use Offset by size of library
  # https://stats.stackexchange.com/questions/66791/where-does-the-offset-go-in-poisson-negative-binomial-regression
  # fit GLM for a row of a sparse matrix, should save some space?

  # pvalue by deviance goodness of fit: https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
  # offset is in log because the model says the log counts is equal to RHS

  if (!is.null(nrow(jrow))){
    # probably a matrix of many rows, sum them up
    print(paste("Merging", nrow(jrow), "rows"))
    row <- Matrix::colSums(jrow)
  }
  dat <- data.frame(cell = cnames, ncuts = jrow, stringsAsFactors = FALSE) %>%
    left_join(., dat.annots.filt.mark, by = "cell") %>%
    left_join(., ncuts.cells.mark, by = "cell")

  m1.pois <- glm(ncuts ~ 1 + ptime + offset(log(ncuts.total)), data = dat, family = "poisson")
  mnull.pois <- glm(ncuts ~ 1 + offset(log(ncuts.total)), data = dat, family = "poisson")


  if (!returnobj){
    jsum <- anova(mnull.pois, m1.pois)
    pval <- pchisq(jsum$Deviance[[2]], df = jsum$Df[[2]], lower.tail = FALSE)

    if (!with.se){
      out.dat <- data.frame(pval = pval,
                            dev.diff = jsum$Deviance[[2]],
                            df.diff = jsum$Df[[2]],
                            t(as.data.frame(coefficients(m1.pois))),
                            stringsAsFactors = FALSE)
    } else {
      estimates <- summary(m1.pois)$coefficients[, "Estimate"]
      names(estimates) <- make.names(paste(names(estimates), ".Estimate", sep = ""))
      stderrors <- summary(m1.pois)$coefficients[, "Std. Error"]
      names(stderrors) <- make.names(paste(names(stderrors), ".StdError", sep = ""))
      out.dat <- data.frame(pval = pval,
                            dev.diff = jsum$Deviance[[2]],
                            df.diff = jsum$Df[[2]],
                            t(as.data.frame(c(estimates, stderrors))),
                            stringsAsFactors = FALSE)
    }

    if (!is.null(jbin)){
      out.dat$bin <- jbin
      rownames(out.dat) <- jbin
    }
    return(out.dat)
  } else {
    return(list(fit.full = m1.pois, fit.null = mnull.pois, dat.input = dat))
  }
}





# Constants ---------------------------------------------------------------

ncores <- 32
outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs"
jprefix <- "var_filtered"
jsuffix <- "manual2nocenter_K36_K9m3_K36-K9m3"
jname <- paste(jprefix, jsuffix, sep = "_")
outdir <- file.path(outmain, jsuffix)
dir.create(outdir)

outf <- file.path(outdir, paste0("glm_poisson_fits_output.", jsuffix, ".", Sys.Date(), ".RData"))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Check whether run on projection or on mixed -----------------------------

# projection
inf.meta.proj <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/celltyping_K9m3_first_try.2021-08-02.txt"
dat.meta.proj <- fread(inf.meta.proj)

ggplot(dat.meta.proj, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# mixed
inf.mixed <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K9m3.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K9m3.K-30.Robj"
load(inf.mixed, v=T)

tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cell2ctype <- hash::hash(dat.meta.proj$cell, dat.meta.proj$cluster)

# compare with dat.umap
dat.umap$ctype <- sapply(dat.umap$cell, function(x) cell2ctype[[x]])

ggplot(dat.umap, aes(x = umap1, y = umap2, color = ctype)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Fit on the unmixed data  -----------------------------------------------

# use princurve
library(princurve)

# fit on UMAP
clsts.keep <- c("Early", "Intermediate1", "Intermediate2", "Late")

dat.umap.filter <- subset(dat.umap, ctype %in% clsts.keep) %>%
  filter(umap2 > -0.5) %>%
  ungroup() %>%
  mutate(umap2 = 0.1 * scale(umap2, center = TRUE, scale = TRUE),
         umap1 = scale(umap1, center = TRUE, scale = TRUE))

mat.filter <- as.matrix(subset(dat.umap.filter, select = c(umap1, umap2)))
rownames(mat.filter) <- dat.umap.filter$cell

pcurveout <- princurve::principal_curve(x = as.matrix(mat.filter))

pcurve.dat <- data.frame(pcurveout$s, cell = rownames(pcurveout$s), ptime = pcurveout$lambda, stringsAsFactors = FALSE)

ggplot() +
  geom_point(mapping = aes(x = umap1, y = umap2), data = dat.umap.filter) +
  geom_point(mapping = aes(x = umap1, y = umap2, color = ptime), data = pcurve.dat) +
  scale_color_viridis_c() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cell2ptime <- hash::hash(pcurve.dat$cell, pcurve.dat$ptime)

dat.umap$ptime <- sapply(dat.umap$cell, function(x) AssignHash(x, cell2ptime, NA))

ggplot(dat.umap, aes(x = umap1, y = umap2, color = ptime)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Fit genes that follow ptime  --------------------------------------------------------------

# load raw counts
cells.keep <- subset(dat.umap, !is.na(ptime))$cell
count.mat.filt <- count.mat[, cells.keep]

# fit poisson regression
dat.annots.filt <- subset(dat.umap, cell %in% cells.keep, select = c(cell, ptime))
dat.annots.filt.othercols <- subset(dat.umap, cell %in% cells.keep, select = c(-ptime))
ncuts.cells <- data.frame(cell = colnames(count.mat.filt), ncuts.total = colSums(count.mat.filt), stringsAsFactors = FALSE)



# jrow.names <- rownames(count.mat.filt)
# names(jrow.names) <- jrow.names
#
# system.time(
#   jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
#     jrow <- jmat.mark[jrow.name, ]
#     jout <- FitGlmRowPtime.withse(jrow = jrow, cnames = cells.keep, dat.annots.filt.mark = dat.annots.filt, ncuts.cells.mark = ncuts.cells, jbin =jbin , returnobj = FALSE, with.se = TRUE)
#     return(jout)
#   }, mc.cores = ncores)
# )
#
# # Ssave outputs -----------------------------------------------------------
#
# save(jfits.lst, dat.annots.filt, ncuts.cells, count.mat.filt, file = outf)

# check output
out.check <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/manual2nocenter_K36_K9m3_K36-K9m3/glm_poisson_fits_output.manual2nocenter_K36_K9m3_K36-K9m3.2021-08-04.RData"
load(out.check, v=T)

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
  geom_hline(yintercept = 0, linetype = "dotted", size = 1, color = "red")

# what are the top hits?

head(print(dat.slope.pval %>% arrange(pval)))

ggplot(dat.annots.filt %>% left_join(., dat.annots.filt.othercols), aes(x = umap1, y = umap2, color = ptime)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.slope.pval.sorted <- dat.slope.pval %>%
  arrange(pval) %>%
  mutate(qval = p.adjust(pval))

print(head(dat.slope.pval.sorted))

jbin <- "chr3:106150000-106200000"
jbin <- "chr9:3000000-3050000"
jbin <- dat.slope.pval.sorted$bin[5]

jbin <- dat.slope.pval.sorted$bin[3]
jbin <- dat.slope.pval.sorted$bin[4]
jbin <- dat.slope.pval.sorted$bin[5]
jbin <- dat.slope.pval.sorted$bin[6]
jbin <- dat.slope.pval.sorted$bin[7]

jbin <- dat.slope.pval.sorted$bin[1]
jbin <- dat.slope.pval.sorted$bin[2]

jcheck <- data.frame(counts = count.mat.filt[jbin, ], cell = colnames(count.mat.filt), stringsAsFactors = FALSE) %>%
  left_join(., dat.annots.filt)  %>%
  left_join(., subset(dat.umap, select = c(cell, umap1, umap2))) %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]])


ggplot(jcheck, aes(x = ptime, y = counts)) +
  geom_point(alpha = 0.1) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot on UMAP

ggplot(jcheck, aes(x = umap1, y = umap2, color = log(counts + 1))) +
  theme_bw() +
  facet_wrap(~stage) +
  geom_point() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

