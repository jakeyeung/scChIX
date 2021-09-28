# Jake Yeung
# Date of Creation: 2021-08-12
# File: ~/projects/scChIX/analysis_scripts/pseudotime/3-fit_celltypes_K36.multicore.discrete.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

ncores <- 32

# Functions ---------------------------------------------------------------


FitGlmRowClusters.withse <- function(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = NULL, returnobj=FALSE, with.se = FALSE){
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

  m1.pois <- glm(ncuts ~ 1 + cluster + offset(log(ncuts.total)), data = dat, family = "poisson")
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



hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "K36"
jname <- "manual2nocenter_K36_K9m3_K36-K9m3"

outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/by_clusters", jname)
dir.create(outdir)
outf <- file.path(outdir, paste0("glm_poisson_fits_output.clusters.", jname, ".", jmark, ".RData"))

# Load rawcounts  --------------------------------------------------------------

inf.obj <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_", jname, "/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K36.K-30.Robj"))
assertthat::assert_that(file.exists(inf.obj))

load(inf.obj, v=T)

# count.mat[1:5, 1:5]

# # marks sometimes contain other marks? doublecheck
# inf.mat <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables_dedup_from_merged/counts_tables_10000/gastru_merged_", jmark, ".rowsdeduped.tagged.sorted.countTable.binsize_10000.csv"))
# mat <- ReadMatSlideWinFormat(inf.mat)

# Load celltype annotations -----------------------------------------------

inf.ctype <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots/metadata_K36-K9me3.", jmark, ".2021-08-09.txt"))
dat.ctype <- fread(inf.ctype)

cells.keep <- dat.ctype$cell

ggplot(dat.ctype, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~type) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# make Epithelial the reference celltype

dat.ctype <- dat.ctype %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "Epithelial", "aEpithelial", cluster),
         cluster = ifelse(startsWith(cluster, "NeuralTubeNeuralProgs"), "NeuralTubeNeuralProgs", cluster))

ggplot(dat.ctype, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~type) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# rename NeuralTubeNeuralProgs2 and 3 to no numbers


# Get norm constants ------------------------------------------------------


# Fit each gene  ---------------------------------------------------------------


cnames <- colnames(count.mat)
dat.annots.filt.mark <- dat.ctype
ncuts.cells.mark <- data.frame(cell = colnames(count.mat), ncuts.total = colSums(count.mat), stringsAsFactors = FALSE)

# test on one
# set.seed(0)
# jrow.i <- sample(seq_len(nrow(count.mat)), size = 1)
# jbin <- rownames(count.mat)[jrow.i]
# jrow <- count.mat[jrow.i, ]
# jfit.out <- FitGlmRowClusters.withse(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jbin, returnobj = FALSE, with.se = TRUE)

# fit all


jrow.names <- rownames(count.mat)
names(jrow.names) <- jrow.names

print("fitting genes")
system.time(
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- count.mat[jrow.name, ]
    # jout <- scChIX::FitGlmRowPtime.withse(jrow = jrow, cnames = cells.keep, dat.annots.filt.mark = dat.annots.filt, ncuts.cells.mark = ncuts.cells, jbin = jrow.name , returnobj = FALSE, with.se = TRUE)
    jout <- FitGlmRowClusters.withse(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jrow.name, returnobj = FALSE, with.se = TRUE)
    return(jout)
  }, mc.cores = ncores)
)

save(jfits.lst, dat.annots.filt.mark, ncuts.cells.mark, count.mat, dat.ctype, file = outf)



