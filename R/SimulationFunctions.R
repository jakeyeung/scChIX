# Jake Yeung
# Date of Creation: 2021-09-20
# File: ~/projects/scChIX/R/SimulationFunctions.R
# Functions for simulating scChIX data


MeanAcrossClusters <- function(count.mat, cnames.keep.lst, jfunc = rowMeans){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    jfunc(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}

MergeMarksEstimateOverlap <- function(bin.dat){
  # bin.dat <- ctype.sim.counts.lst$A$bin.data
  bin.dat <- bin.dat[gtools::mixedorder(bin.dat$Bin), ]
  # mat.prob <- mat.prob[, grepl(grepsuffix, colnames(mat.prob))]
  # mat.dbl <- mats$`mark1-mark2`[, grepl(grepsuffix, colnames(mats$`mark1-mark2`))]

  # compare with ground truth
  jnbins <- nrow(bin.dat)
  frac.swap <- ctype.params.lst$A$frac.mutual.excl / 2
  bottom.bins <- (bin.dat %>% dplyr::arrange(BinMean))$Bin[1:(jnbins * frac.swap)]
  top.bins <- (bin.dat %>% dplyr::arrange(desc(BinMean)))$Bin[1:(jnbins * frac.swap)]


  bin.means.A1 <- data.frame(mark = "mark1", BinMean = bin.dat$BinMean, annot = bin.dat$annot, stringsAsFactors = FALSE)
  rownames(bin.means.A1) <- bin.dat$Bin

  bin.means.A2 <- data.frame(mark = "mark2", BinMean = bin.dat$BinMean, stringsAsFactors = FALSE)
  rownames(bin.means.A2) <- bin.dat$Bin

  bin.means.A2 <- scChIX::SwapBins(ctype1.sim.counts.b = bin.means.A2, top.bins = top.bins, bottom.bins = bottom.bins, frac.swap = frac.swap)

  bin.means.A1$Bin <- rownames(bin.means.A1)
  bin.means.A2$Bin <- rownames(bin.means.A2)

  bin.means.A.merged <- left_join(bin.means.A1, bin.means.A2, by = "Bin") %>%
    rowwise() %>%
    mutate(BinMean.dbl = (BinMean.x) / (BinMean.x + BinMean.y))
  return(bin.means.A.merged)
}



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

  print(Matrix::nnzero(ctype1.sim.counts) / length(ctype1.sim.counts))


  bin.data <- as.data.frame(rowData(ctype1.sim))
  bin.data$Bin <- as.character(bin.data$Bin)
  bin.data$Bin.Orig <- as.character(bin.data$Bin)
  # shuffle bins
  bin.data$Bin <- sample(bin.data$Bin.Orig, size = nrow(bin.data), replace = FALSE)

  assertthat::assert_that(all(sort(unique(bin.data$Bin)) == sort(unique(bin.data$Bin.Orig))))

  print(head(bin.data %>% dplyr::arrange(desc(BinMean))))

  print(str(ctype1.sim.counts))

  # bin.counts <- sort(rowSums(ctype1.sim.counts), decreasing = TRUE)
  # print(head(bin.counts))

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

