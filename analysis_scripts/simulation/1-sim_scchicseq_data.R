# Jake Yeung
# Date of Creation: 2021-09-17
# File: ~/projects/scChIX/analysis_scripts/simulation/sim_scchicseq_data.R
#

rm(list=ls())


# library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(irlba)
library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

# library(devtools)
# install_github("bowang-lab/simATAC")

library(simATAC)


# Functions ---------------------------------------------------------------


ConcatMat <- function(ctype.sim.counts.lst, ctypes, jmark){
  mat.mark1.lst <- lapply(ctypes, function(jctype){
    ctype.sim.counts.lst[[jctype]][[jmark]]
  })
  rnames.common <- unique(gtools::mixedsort(unlist(lapply(mat.mark1.lst, rownames))))

  mat.mark1.lst <- lapply(mat.mark1.lst, function(jmat){
    mat.rearranged <- jmat[rnames.common, ]
  })
  mat.mark1 <- do.call(cbind, mat.mark1.lst)
  return(mat.mark1)
}



SwapBins <- function(ctype1.sim.counts.b, top.bins, bottom.bins, frac.swap, jnbins){
  ctype1.sim.counts.b.swapped <- ctype1.sim.counts.b
  ctype1.sim.counts.b.swapped[bottom.bins, ] <- ctype1.sim.counts.b[top.bins, ]
  ctype1.sim.counts.b.swapped[top.bins, ] <- ctype1.sim.counts.b[bottom.bins, ]
  return(ctype1.sim.counts.b.swapped)
}

SimulateChICseq <- function(ctype.name = "A", jseed = 0, shuffle.rows.seed = 0, frac.mutual.excl = 0.5, jspec = "mm10",
                            jncells.per.rep = 250, jnbins = 5000, jlibmean = 12,
                            jlibsd = 1, jlibp = 0.5, jzp = 1){


  nreps <- 3  # two singles, one doubleincubation mark
  frac.swap <- frac.mutual.excl / 2  # half mutual exclusive bins aree present in A, other half present in B

  ctype1 <- newsimATACCount()
  ctype1 <- setParameters(ctype1, species = jspec)
  ctype1 <- setParameters(ctype1, seed = jseed)
  ctype1 <- setParameters(ctype1, nCells = jncells.per.rep * nreps)
  ctype1 <- setParameters(ctype1, nBins = jnbins)
  ctype1 <- setParameters(ctype1, lib.mean1 = jlibmean)
  ctype1 <- setParameters(ctype1, lib.mean2 = jlibmean)
  ctype1 <- setParameters(ctype1, lib.sd1 = jlibsd)
  ctype1 <- setParameters(ctype1, lib.sd1 = jlibsd)
  ctype1 <- setParameters(ctype1, lib.prob = jlibp)
  ctype1 <- setParameters(ctype1, non.zero.pro = jzp)


  ctype1.sim <- simATACSimulate(ctype1)

  plot(density(log10(ctype1.sim$LibSize)))

  ctype1.sim.counts <- counts(ctype1.sim)

  Matrix::nnzero(ctype1.sim.counts) / length(ctype1.sim.counts)


  bin.data <- as.data.frame(rowData(ctype1.sim))
  bin.data$Bin <- as.character(bin.data$Bin)
  bin.data$Bin.Orig <- as.character(bin.data$Bin)
  # shuffle bins
  bin.data$Bin <- sample(bin.data$Bin.Orig, size = nrow(bin.data), replace = FALSE)

  assertthat::assert_that(all(sort(unique(bin.data$Bin)) == sort(unique(bin.data$Bin.Orig))))

  print(head(bin.data %>% dplyr::arrange(desc(BinMean))))

  bin.counts <- sort(rowSums(ctype1.sim.counts), decreasing = TRUE)
  print(head(bin.counts))

  # shuffle bins
  ctype1.sim.counts.orig <- ctype1.sim.counts

  # shuffle rows to create cell type-specific matrices
  rownames(ctype1.sim.counts) <- bin.data$Bin

  colnames(ctype1.sim.counts) <- paste(colnames(ctype1.sim.counts), ctype.name, sep = "_")
  print(head(ctype1.sim.counts[1:5, 1:5]))


  # Annotate bins -----------------------------------------------------------

  bottom.bins <- (bin.data %>% dplyr::arrange(BinMean))$Bin[1:(jnbins * frac.swap)]
  top.bins <- (bin.data %>% dplyr::arrange(desc(BinMean)))$Bin[1:(jnbins * frac.swap)]
  mutual.excl.bins <- c(bottom.bins, top.bins)
  common.bins <- bin.data$Bin[!bin.data$Bin %in% mutual.excl.bins]

  bin.data$annot <- sapply(bin.data$Bin, function(b) ifelse(b %in% common.bins, "Common", "MutualExcl"))

  # Split counts into two  --------------------------------------------------

  ctype1.sim.counts.a <- ctype1.sim.counts[, 1:jncells.per.rep]
  ctype1.sim.counts.dbl.a <- ctype1.sim.counts[, (jncells.per.rep * 2 + 1) : (jncells.per.rep * 3)]
  # these mats need bin swapping to map to "mark b"
  ctype1.sim.counts.b <- ctype1.sim.counts[, (jncells.per.rep + 1) : (jncells.per.rep * 2)]
  ctype1.sim.counts.dbl.b <- ctype1.sim.counts.dbl.a   # same as dbl.a because technical effects are controlled within a cell

  # swap top X percent with bottom X percent from (a) to (b)
  ctype1.sim.counts.b.swapped <- SwapBins(ctype1.sim.counts.b, top.bins = top.bins, bottom.bins = bottom.bins, frac.swap = frac.swap, jnbins = jnbins)
  assertthat::assert_that(all(rownames(ctype1.sim.counts.b.swapped) == rownames(ctype1.sim.counts.a)))
  # plot(rowSums(ctype1.sim.counts.a[mutual.excl.bins, ]), rowSums(ctype1.sim.counts.b.swapped[mutual.excl.bins, ]))
  # plot(rowSums(ctype1.sim.counts.a[common.bins, ]), rowSums(ctype1.sim.counts.b.swapped[common.bins, ]))


  # # these bins are obviously correlated
  # plot(rowSums(ctype1.sim.counts.a), rowSums(ctype1.sim.counts.b))


  # swap top X percent with bottom X percent from (a) to (b)
  ctype1.sim.counts.dbl.b.swapped <- SwapBins(ctype1.sim.counts.dbl.b, top.bins = top.bins, bottom.bins = bottom.bins, frac.swap = frac.swap, jnbins = jnbins)
  assertthat::assert_that(all(rownames(ctype1.sim.counts.dbl.b.swapped) == rownames(ctype1.sim.counts.dbl.a)))

  ctype1.sim.counts.dbl.sum <- ctype1.sim.counts.dbl.a + ctype1.sim.counts.dbl.b.swapped
  colnames(ctype1.sim.counts.dbl.sum) <- paste(colnames(ctype1.sim.counts.dbl.a),
                                           colnames(ctype1.sim.counts.dbl.b.swapped), sep = "x")

  ctype1.sim.counts.lst <- list(mark1 = ctype1.sim.counts.a, mark2 = ctype1.sim.counts.b.swapped, markdbl = ctype1.sim.counts.dbl.sum, bin.data = bin.data)
  return(ctype1.sim.counts.lst)
}


# Constants ---------------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_inputs/countmats")

jspec <- "mm10"
jseed <- 0
jncells.per.rep <- 250
jnbins <- 5000
jlibmean <- 12
jlibsd <- 1
jlibp <- 0.5
jzp <- 1
# frac.mutual.excl <- 0.5
shuffle.rows.seed <- jseed

sim.params <- list(jspec = jspec,
                   jseed = jseed,
                   jncells.per.rep = jncells.per.rep,
                   jnbins = jnbins,
                   jlibmean = jlibmean,
                   jlibsd = jlibsd,
                   jlibp = jlibp,
                   jzp = jzp)


# Run simulation ----------------------------------------------------------

ctypes <- c("A", "B", "C")
names(ctypes) <- ctypes

ctype.params <- list(ctype.name = c("A", "B", "C"),
                     jseed = c(0, 123, 999),
                     shuffle.rows.seed = c(0, 123, 999),
                     frac.mutual.excl = c(0.5, 0.5, 0.5))

ctype.params.lst <- lapply(ctypes, function(ctype){
  i <- which(ctype.params$ctype.name == ctype)
  return(list(ctype = ctype,
              jseed = ctype.params$jseed[[i]],
              shuffle.rows.seed = ctype.params$shuffle.rows.seed[[i]],
              frac.mutual.excl = ctype.params$frac.mutual.excl[[i]]))
})

ctype.sim.counts.lst <- lapply(ctypes, function(jctype){
  SimulateChICseq(ctype.name = jctype,
                  jseed = ctype.params.lst[[jctype]]$jseed,
                  shuffle.rows.seed = ctype.params.lst[[jctype]]$shuffle.rows.seed,
                  frac.mutual.excl = ctype.params.lst[[jctype]]$frac.mutual.excl)
})


# Create simulated matrix -------------------------------------------------


jmark1 <- "mark1"
jmark2 <- "mark2"
jmarkdbl <- paste(jmark1, jmark2, sep = "-")

mat.mark1 <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "mark1")
mat.mark2 <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "mark2")
mat.markdbl <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "markdbl")

mat.mark.lst <- list(mat.mark1, mat.mark2, mat.markdbl)
names(mat.mark.lst) <- c(jmark1, jmark2, jmarkdbl)

jmarks <- names(mat.mark.lst)
names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 150
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

dat.umap.lst <- lapply(mat.mark.lst, function(jmat){
  print(head(jmat[1:5, 1:5]))
  lsi.out <- scchicFuncs::RunLSI(as.matrix(jmat))
  dat.umap.mark <- scchicFuncs::DoUmapAndLouvain(topics.mat = lsi.out$u, jsettings = jsettings)
  return(dat.umap.mark)
})

m.lst <- lapply(dat.umap.lst, function(jdat){
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

print(m.lst)

# Write to output ---------------------------------------------------------

sim.files <- file.path(outdir, paste0("ATAC_simulator_params_and_outputs.RData"))

save(ctype.params.lst, ctype.sim.counts.lst, file = sim.files)

# save each mark
for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("countmat_var_filt.", jmark, ".rds"))
  saveRDS(mat.mark.lst[[jmark]], file = outf)
}



# Checks  -----------------------------------------------------------------


dat.meta.lst <- lapply(jmarks, function(jmark){
  dat.meta.tmp <- data.frame(cell = colnames(mat.mark.lst[[jmark]]), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(ctype = substr(cell, start = nchar(cell), nchar(cell)))
})

cnames.keep.lst.lst <- lapply(jmarks, function(jmark){
  jlist <- split(x = dat.meta.lst[[jmark]], f = dat.meta.lst[[jmark]]$ctype)
  lapply(jlist, function(x) x$cell)
})

mats.pbulk <- lapply(jmarks, function(jmark){
  vec.lst <- SumAcrossClusters(mat.mark.lst[[jmark]], cnames.keep.lst = cnames.keep.lst.lst[[jmark]])
  dat <- as.data.frame(vec.lst)
})


plot(mats.pbulk$mark1$A, mats.pbulk$mark2$A)

plot(mats.pbulk$mark1$A, mats.pbulk$mark1$C)

plot(mats.pbulk$mark2$A, mats.pbulk$mark2$B)
plot(mats.pbulk$mark2$A, mats.pbulk$mark2$C)

plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$B)
plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$C)



plot(mats.pbulk$mark1$A + mats.pbulk$mark2$A, mats.pbulk$`mark1-mark2`$A)

plot(mats.pbulk$mark1$B + mats.pbulk$mark2$B, mats.pbulk$`mark1-mark2`$B)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C, log = "xy")

