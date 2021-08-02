# Jake Yeung
# Date of Creation: 2021-07-10
# File: ~/projects/scChIX/analysis_scripts/3e-get_top_bins_and_genes.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load bins and TSSs ------------------------------------------------------

jmarks <- c("K36", "K9m3")
names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# load the two ldas

inf.lda.k36 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_TES_k36_cleaned/lda_outputs.countmat_TES_cleaned.K36.2021-07-09.K-30.binarize.FALSE/ldaOut.countmat_TES_cleaned.K36.2021-07-09.K-30.Robj")
inf.lda.k9 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_dbl_cleaned/lda_outputs.countmat_output_filt.K9m3.2021-07-09.K-30.binarize.FALSE/ldaOut.countmat_output_filt.K9m3.2021-07-09.K-30.Robj")

load(inf.lda.k36, v=T)
out.lda.k36 <- out.lda
count.mat.k36 <- count.mat
load(inf.lda.k9, v=T)
out.lda.k9 <- out.lda
count.mat.k9 <- count.mat


# Get top 500 bins for every topic  ---------------------------------------

tm.result.k36 <- posterior(out.lda.k36)
tm.result.k9 <- posterior(out.lda.k9)

dat.umap.k36 <- DoUmapAndLouvain(topics.mat = tm.result.k36$topics, jsettings = jsettings)
dat.umap.k9 <- DoUmapAndLouvain(topics.mat = tm.result.k9$topics, jsettings = jsettings)

m.k36 <- ggplot(dat.umap.k36, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.k9 <- ggplot(dat.umap.k9, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.k36)
print(m.k9)


# Get top bins  -----------------------------------------------------------

topbins.keep <- 250

GetTopFeatureNames <- function(jrow, cnames, jtopbins.keep){
  names(jrow) <- cnames
  jrow.rank <- rank(-1 * jrow)
  jrow.rank.keep <- jrow.rank < topbins.keep
  features.keep <- names(jrow.rank.keep)[jrow.rank.keep]
  return(features.keep)
}

tm.result.lst <- list("K36" = tm.result.k36,
                      "K9m3" = tm.result.k9)

top.features.common.lst <- lapply(tm.result.lst, function(tm.result){
  cnames <- colnames(tm.result$terms)
  top.features.lst <- apply(tm.result$terms, 1, function(jrow) GetTopFeatureNames(jrow, cnames, topbins.keep))
  top.features.common <- unique(unlist(top.features.lst))
  return(top.features.common)
})



# Write new count mat, filtered by features  ------------------------------

count.mat.lst <- list("K36" = count.mat.k36,
                      "K9m3" = count.mat.k9)

count.mat.filt.lst <- lapply(jmarks, function(jmark){
  features.keep <- top.features.common.lst[[jmark]]
  count.mat.tmp <- count.mat.lst[[jmark]][features.keep, ]
})


# Write output  -----------------------------------------------------------

# count.mats.filt and bedfiles of all the locations
coords.k36 <- sapply(top.features.common.lst$K36, function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)
coords.k9 <- sapply(top.features.common.lst$K9m3, function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)
coords.k9 <- gsub("^chr", "", coords.k9)

coords.merged <- c(coords.k36, coords.k9)
chromos <- sapply(coords.merged, function(x) GetChromo(x), USE.NAMES = FALSE)
starts <- sapply(coords.merged, function(x) GetStart(x), USE.NAMES = FALSE)
ends <- sapply(coords.merged, function(x) GetEnd(x), USE.NAMES = FALSE)
namesvec <- c(top.features.common.lst$K36, top.features.common.lst$K9m3)

dat.coords <- data.frame(Chromo = chromos, Start = starts, End = ends, Name = namesvec, stringsAsFactors = FALSE)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/coords_filtered"
outbed <- file.path(outdir, paste0("coords_K36_genebodies_K9m3_bins.", Sys.Date(), ".bed"))
fwrite(dat.coords, file = outbed, sep = "\t", col.names = FALSE)

outmat1 <- file.path(outdir, paste0("countmat_featurefilt.", jmarks[[1]], ".", Sys.Date(), ".rds"))
outmat2 <- file.path(outdir, paste0("countmat_featurefilt.", jmarks[[2]], ".", Sys.Date(), ".rds"))

# write mats
saveRDS(count.mat.filt.lst[[jmarks[[1]]]], file = outmat1)
saveRDS(count.mat.filt.lst[[jmarks[[2]]]], file = outmat2)

