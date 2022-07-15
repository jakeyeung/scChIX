# Jake Yeung
# Date of Creation: 2021-09-20
# File: ~/projects/scChIX/R/SimulationFunctions.R
# Functions for simulating scChIX data


#' Wrangling function: Take count mat and summarize into mean across selected columns (eg by columns that are from same cluster)
#' 
#' @param count.mat Count matrix, usually rows are genes and columns are cells
#' @param cnames.keep.lst List where names are cluster identities, values in list are column names that are to be summarized
#' @param jfunc Summarizing function, default is rowMeans
#' @return List of count vectors with same names as cnames.keep.lst. Count values are across rows of the matrix 
#' @examples
#' cnames.keep.lst <- lapply(split(dat.final.annots, f = dat.final.annots$cluster), function(x) x$cell) # create cluster to cell input list
#' prob.mat.byclst.lst <- MeanAcrossClusters(prob.mat, cnames.keep.lst = cnames.keep.lst) # average probabilities across clusters
#' @export
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


#' Wrangling function: Takes different cell types from a histone mark X SimulateChICseq and concatenates them together
#' 
#' @param ctype.sim.counts.lst Output from SimulateChICseq which contains different cell types across different histone mark conditions. 
#' @param ctypes Vector of name of cell types to be concatenated together. 
#' @param jmark Name of histone mark from which cell types will be extracted.
#' @return mat.mark1 Matrix of counts, rows are bins, columns are cells.
#' @examples
#' ctype.sim.counts.lst <- lapply(ctypes, function(jctype){
#' SimulateChICseq(ctype.name = jctype,
#'                 jseed = ctype.params.lst[[jctype]]$jseed,
#'                 shuffle.rows.seed = ctype.params.lst[[jctype]]$shuffle.rows.seed,
#'                 frac.mutual.excl = ctype.params.lst[[jctype]]$frac.mutual.excl)
#' })
#' mat.mark1 <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "mark1")
#' @export
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

#' Simulate count matrices from a Poisson distribution
#' 
#' @param ctype.name Name of "celltype" that these single counts come from. Default "A"
#' @param jseed Random seed for simulating count matrices from Poisson distribution
#' @param shuffle.rows.seed Random seed for shuffling rows. Changing this number and rerunning would generate counts from another "celltype".
#' @param frac.mutual.excl Fraction of bins that will be made mutually exclusive between the second mark. This is done by swapping bins with the highest counts (i.e. largest lambda from Poisson) with the lowest counts.
#' @param jspec Species is mm10 or hg38. Defines position information.
#' @param jncells.per.rep number of cells per replicate. Generates 3 replicates. One for histone mark 1, one for histone mark 2, and one with histone mark 1 + 2
#' @param jnbins Number of bins in the data
#' @param jlibmean Average (mean) library size
#' @param libsd Standard deviation in library size
#' @param jzp Additional dropouts from Poisson (will still be Poisson after dropouts)
#' @return ctype1.sim.counts.lst Simulated counts from the three replicates (hist 1, hist 2, hist1+2) and annotations of whether the bin is overlapping or mutually exclusive.
#' @examples
#' 
#' jspec <- "hg38"
#' jseed <- 0
#' jncells.per.rep <- 250
#' jnbins <- 10000
#' jlibmean <- 12
#' jlibsd <- 1
#' jlibp <- 0.5
#' jzp <- 1
#' shuffle.rows.seed <- jseed
#' sim.params <- list(jspec = jspec,
#' jseed = jseed,
#' jncells.per.rep = jncells.per.rep,
#' jnbins = jnbins,
#' jlibmean = jlibmean,
#' jlibsd = jlibsd,
#' jlibp = jlibp,
#' jzp = jzp)
#' ctype.sim.counts.lst <- lapply(ctypes, function(jctype){
#' ctypes <- c("A", "B", "C")
#' names(ctypes) <- ctypes
#' 
#' ctype.params <- list(ctype.name = c("A", "B", "C"),
#'                      jseed = c(0, 123, 999),
#'                      shuffle.rows.seed = c(0, 123, 999),
#'                      frac.mutual.excl = c(0.5, 0.5, 0.5))
#' 
#' ctype.params.lst <- lapply(ctypes, function(ctype){
#'   i <- which(ctype.params$ctype.name == ctype)
#'   return(list(ctype = ctype,
#'               jseed = ctype.params$jseed[[i]],
#'               shuffle.rows.seed = ctype.params$shuffle.rows.seed[[i]],
#'               frac.mutual.excl = ctype.params$frac.mutual.excl[[i]]))
#' })
#' SimulateChICseq(ctype.name = jctype,
#'                 jseed = ctype.params.lst[[jctype]]$jseed,
#'                 shuffle.rows.seed = ctype.params.lst[[jctype]]$shuffle.rows.seed,
#'                 frac.mutual.excl = ctype.params.lst[[jctype]]$frac.mutual.excl)
#' })
#' @export
SimulateChICseq <- function(ctype.name = "A", jseed = 0, shuffle.rows.seed = 0, frac.mutual.excl = 0.5, jspec = "hg38",
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
  set.seed(shuffle.rows.seed)
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

