# Jake Yeung
# Date of Creation: 2021-09-20
# File: ~/projects/scChIX/analysis_scripts/simulation/4-plot_technical_factors_of_simulation_input_data.R
# Find convincing plots to show the input simulation matrix makes sense

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scChIX)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

library(GGally)



# Change this  ------------------------------------------------------------

frac.mutexcl <- 0.5
frac.mutexcl.str <- sprintf("%.1f", round(frac.mutexcl, 1))

# frac.mutexcl <- 0.01
# frac.mutexcl.str <- sprintf("%.2f", round(frac.mutexcl, 2))


# Constants ---------------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("mark1", "mark2", "mark1-mark2"); names(jmarks) <- jmarks
ctypes <- c("A", "B", "C"); names(ctypes) <- ctypes

inmain <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data")
indir <- file.path(inmain, paste("frac_mutexcl", frac.mutexcl.str, sep = "_"))
assertthat::assert_that(dir.exists(indir))

outpdf <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/simulation_data_downstream/plots/technical_plots.FracMutExcl_", frac.mutexcl.str, ".", Sys.Date(), ".pdf"))


# Load sim outputs -------------------------------------------------------

inf.sim <- file.path(indir, "snakemake_inputs/countmats/ATAC_simulator_params_and_outputs.RData")
load(inf.sim, v=T)



# Load raw data  ----------------------------------------------------------


inf.mats <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, "snakemake_inputs", "countmats", paste0("countmat_var_filt.", jmark, ".rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

mats.lst <- lapply(inf.mats, function(jinf){
  readRDS(jinf)
})

dat.meta.lst <- lapply(jmarks, function(jmark){
  jmat <- mats.lst[[jmark]]
  ctypes.vec <- sapply(colnames(jmat), function(x) substr(x, start = nchar(x), stop = nchar(x)))
  dat.meta.tmp <- data.frame(cell = colnames(jmat), ctype = ctypes.vec, libsize = colSums(jmat), stringsAsFactors = FALSE)
  cells.sparsity <- apply(jmat, MARGIN = 2, FUN = function(jcol){
    1 - Matrix::nnzero(jcol) / length(jcol)
  })
  dat.sparsity <- data.frame(cell = names(cells.sparsity), frac.zeros = cells.sparsity, stringsAsFactors = FALSE)
  dat.meta.out <- left_join(dat.meta.tmp, dat.sparsity, by = "cell")
  dat.meta.out$mark <- jmark
  return(dat.meta.out)
})


# Calculate technical factors ---------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# library
m.libsize <- ggplot(dat.meta.lst %>% bind_rows(), aes(x = libsize, fill = mark)) +
  geom_density(alpha = 0.25) +
  scale_x_log10() +
  xlab("Library size per cell") +
  scale_fill_manual(values = cbPalette) +
  theme_bw() +
  ggtitle("Library size across cells") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# sparsity
m.sparse <- ggplot(dat.meta.lst %>% bind_rows(), aes(x = frac.zeros, fill = mark)) +
  geom_density(alpha = 0.25) +
  xlab("Fraction of Empty Bins") +
  scale_fill_manual(values = cbPalette) +
  theme_bw() +
  ggtitle("Sparsity across cells") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# pseudobulk correlation between marks

dat.meta.lst.byctype <- lapply(dat.meta.lst, function(jdat){
  lapply(split(jdat, f = jdat$ctype), function(x) x$cell)
})

mats.pbulk.lst <- lapply(jmarks, function(jmark){
  mat.pbulk <- as.data.frame(SumAcrossClusters(count.mat = mats.lst[[jmark]], cnames.keep.lst = dat.meta.lst.byctype[[jmark]]))
})




# Calculate overlaps  -----------------------------------------------------

bin.means.merged.lst <- lapply(ctypes, function(jctype){
  grepsuffix <- paste0("_", jctype, "$")
  bin.means.merged <- MergeMarksEstimateOverlap(bin.dat = ctype.sim.counts.lst[[jctype]]$bin.data)
  bin.means.merged$ctype <- jctype
  return(bin.means.merged)
})


# m.overlaps <- lapply(ctypes, function(jctype){
#   jdat <- bin.means.merged.lst[[jctype]]
#   ggplot(jdat, aes(x = BinMean.dbl, fill = annot, y =..count../sum(..count..))) +
#     geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
#     theme_bw() +
#     xlab("Fraction of signal belonging to mark1") +
#     ylab("Fraction of bins") +
#     ggtitle(paste(jctype)) +
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# })

m.overlaps <- ggplot(bin.means.merged.lst %>% bind_rows(), aes(x = BinMean.dbl, fill = annot, y =..count../sum(..count..))) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    theme_bw() +
    xlab("Fraction of signal belonging to mark1") +
    ylab("Fraction of bins") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



# Get UMAPs  --------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 150
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

inf.lda.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, paste0("snakemake_outputs/LDA_outputs_init/ldaOut.countmat_var_filt.", jmark, ".Robj"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

out.lda.lst <- lapply(inf.lda.lst, function(jinf){
  load(jinf, v=T)
  return(out.lda)
})

dat.umap.lst <- lapply(out.lda.lst, function(jout){
  DoUmapAndLouvain(posterior(jout)$topics, jsettings = jsettings)
})

dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  dat.umap.lst[[jmark]] %>%
    left_join(., dat.meta.lst[[jmark]])
})

m.umap.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.umap.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})



# Write outputs -----------------------------------------------------------

pdf(outpdf, useDingbats = FALSE)

print(m.libsize)
print(m.sparse)
print(m.overlaps)

print(m.libsize + facet_wrap(~ctype))
print(m.sparse + facet_wrap(~ctype))
print(m.overlaps + facet_wrap(~ctype))

JFuncs::multiplot(m.umap.lst[[1]], m.umap.lst[[2]], m.umap.lst[[3]], cols = 3)

ggpairs(mats.pbulk.lst$mark1)
ggpairs(mats.pbulk.lst$mark2)
ggpairs(mats.pbulk.lst$`mark1-mark2`)

for (ctype in ctypes){
  plot(mats.pbulk.lst$mark1[[ctype]], mats.pbulk.lst$mark2[[ctype]], main = ctype, log = "xy")
}

dev.off()

