# Jake Yeung
# scChIX_pseudotime.R
# 2022-03-11
# DESCRIPTION

FitMixingWeightPtimeLinear.freew <- function (cell.count.raw.merged, jfits.lowess1, jfits.lowess2, 
                                             w.init, w.lower, w.upper, p1.init, p1.lower, p1.upper, p2.init, 
                                             p2.lower, p2.upper, jmethod = "L-BFGS-B", jhessian = FALSE, 
                                             jcontrol = list(ndeps = c(0.01, 0.05, 0.05)))
{
  optim(par = c(w.init, p1.init, p2.init), fn = GetLLMergedPtimeLowessLinear, 
        cell.count.raw.merged = cell.count.raw.merged, fits.lowess1 = jfits.lowess1, 
        fits.lowess2 = jfits.lowess2, method = jmethod, hessian = jhessian, 
        lower = c(w.lower, p1.lower, p2.lower), upper = c(w.upper, 
                                                          p1.upper, p2.upper), 
        control = jcontrol)
}



PredictSignalGenomeWide <- function(loess.fits, ptime, return.y.only = FALSE){
  jgenes.vec <- names(loess.fits); names(jgenes.vec) <- jgenes.vec
  preds <- lapply(jgenes.vec, function(jgene){
    predict(loess.fits[[jgene]], ptime, se = FALSE)
  })
  # create vector
  if (!return.y.only){
    return(unlist(preds))
  } else {
    return(unlist(lapply(preds, function(p) p$y)))
  }
}

PredictSignalGenomeWideLowess <- function(lowess.fits, ptime, linearize = TRUE, normalize = TRUE){
  # assumes log2
  jgenes.vec <- names(lowess.fits); names(jgenes.vec) <- jgenes.vec
  preds <- lapply(jgenes.vec, function(jgene){
    predict(lowess.fits[[jgene]], ptime, se = FALSE)$y
  }) %>%
    unlist()
  # renormalizeA
  if (linearize){
    preds <- 2^preds
  }
  if (normalize){
    preds <- preds / sum(preds)
  }
  return(preds)
}

#' Estimate the signal given a lowess fit and pseudotime.
#' 
#' @param lowess.fits lowess fits for every gene along pseudotime
#' @param ptime pseudotime to estimate the signal at
#' @return preds: for every gene, what is the predicted signal at pseudotime ptime
#' @examples
#' x.raw.unmixed <- parallel::mclapply(all.cells, function(jcell){
#' x.raw <- all.x.raw[[jcell]]
#' ptime1 <- dat.fits.clean.lst[[jcell]]$ptime1
#' ptime2 <- dat.fits.clean.lst[[jcell]]$ptime2
#' p1.cell <- PredictSignalGenomeWideLowessLinear(lowess.fits = lowess.fits.k4me1, ptime = ptime1)
#' p2.cell <- PredictSignalGenomeWideLowessLinear(lowess.fits = lowess.fits.k36me3, ptime = ptime2)
#' x.unmixed.lst <- UnmixRawCounts(x.raw = x.raw, mixweight = w.fixed, p.active = p1.cell, p.repress = p2.cell, random.seed = 0)
#' return(x.unmixed.lst)
#' }, mc.cores = ncores)
#' @export
PredictSignalGenomeWideLowessLinear <- function(lowess.fits, ptime){
  # assumes log2
  jgenes.vec <- names(lowess.fits); names(jgenes.vec) <- jgenes.vec
  preds <- lapply(jgenes.vec, function(jgene){
    predict(lowess.fits[[jgene]], ptime, se = FALSE)$y
  }) %>%
    unlist()
  return(preds)
}


GetLLMergedPtimeLowessLinear <- function(params, cell.count.raw.merged, fits.lowess1, fits.lowess2){
  w <- params[[1]]
  jptime1 <- params[[2]]
  jptime2 <- params[[3]]
  # return neg LL
  w.vec <- c(w, 1-w)  # active proportion, repress proportion
  if (any(w.vec < 0)){
    print(paste("Error message: w is negative:"))
    print(w.vec)
    w.vec[which.min(w.vec)] <- 0
    w.vec[which.max(w.vec)] <- 1
    print("Fixed w:")
    print(w.vec)
  }
  
  prob.vec1 <- PredictSignalGenomeWideLowessLinear(lowess.fits = fits.lowess1, ptime = jptime1)
  prob.vec2 <- PredictSignalGenomeWideLowessLinear(lowess.fits = fits.lowess2, ptime = jptime2)
  merged.prob.vec <- matrixStats::colWeightedMeans(rbind(prob.vec1, prob.vec2), w = w.vec)  # standard weighted.mean is SUPER slow
  
  # print(plot(density(merged.prob.vec)))
  
  # print(range(merged.prob.vec))
  LL <- dmultinom(x = cell.count.raw.merged, prob = merged.prob.vec, log = TRUE)
  # print(sum(merged.prob.vec))
  # print(params)
  # print(LL)
  return(-1 * LL)
}

GetLLMergedPtimeLowessLinear.fixedw <- function(params, cell.count.raw.merged, fits.lowess1, fits.lowess2, w = 0.3){
  # w <- params[[1]]
  jptime1 <- params[[1]]
  jptime2 <- params[[2]]
  # return neg LL
  w.vec <- c(w, 1-w)  # active proportion, repress proportion
  if (any(w.vec < 0)){
    print(paste("Error message: w is negative:"))
    print(w.vec)
    w.vec[which.min(w.vec)] <- 0
    w.vec[which.max(w.vec)] <- 1
    print("Fixed w:")
    print(w.vec)
  }
  
  prob.vec1 <- PredictSignalGenomeWideLowessLinear(lowess.fits = fits.lowess1, ptime = jptime1)
  prob.vec2 <- PredictSignalGenomeWideLowessLinear(lowess.fits = fits.lowess2, ptime = jptime2)
  merged.prob.vec <- matrixStats::colWeightedMeans(rbind(prob.vec1, prob.vec2), w = w.vec)  # standard weighted.mean is SUPER slow
  
  # print(plot(density(merged.prob.vec)))
  
  # print(range(merged.prob.vec))
  LL <- dmultinom(x = cell.count.raw.merged, prob = merged.prob.vec, log = TRUE)
  # print(sum(merged.prob.vec))
  # print(params)
  # print(LL)
  return(-1 * LL)
}


FitMixingWeightPtimeLinear <- function(cell.count.raw.merged, jfits.lowess1, jfits.lowess2, w.init, w.lower, w.upper, p1.init, p1.lower, p1.upper, p2.init, p2.lower, p2.upper, jmethod = "L-BFGS-B", jhessian = FALSE){
  # https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/hessian.nb.html
  # method by Brent is recommended?
  # jmethod <- ""L-BFGS-B""
  # w, ptime, cell.count.raw.merged, fits.lowess1, fits.lowess2
  optim(par = c(w.init, p1.init, p2.init),
        fn = GetLLMergedPtimeLowessLinear,
        cell.count.raw.merged = cell.count.raw.merged,
        fits.lowess1 = jfits.lowess1,
        fits.lowess2 = jfits.lowess2,
        method = jmethod, hessian=jhessian, lower = c(w.lower, p1.lower, p2.lower), upper = c(w.upper, p1.upper, p2.upper), control = list(ndeps = c(0.02, 0.02, 0.02)))
}



#' Infer pseudotime1 and pseudotime2 from a double-incubated cell
#' 
#' @param cell.count.raw.merged vector of raw cell counts from a double-incubated cell. If KNN smoothing, then this vector is the sum of the double-incubated cell and its nearest K neighbors (K=25 used in paper). 
#' @param jfits.lowess1 Named list of "smooth.spline" outputs (lowess fits) for mark1. Each element in list is accessible by a gene name that matches the names in cell.count.raw.merged
#' @param jfits.lowess2 Named list of "smooth.spline" outputs (lowess fits) mark2. Each element in list is accessible by a gene name that matches the names in cell.count.raw.merged
#' @param w.fixed Fixed weighting between mark1 and mark2. In paper we used w=0.77, optimized by running scChIX for different w's. Fixing w for all cells makes for more stable outputs than a cell-specific w
#' @param p1.init init pseudotime1 [0, 1]
#' @param p2.init init pseudotime2 [0, 1]
#' @param p1.lower lowerbound for pseudotime1 (eg 0)
#' @param p2.lower lowerbound for pseudotime2 (eg 0)
#' @param p1.upper upperbound for pseudotime1 (eg 1)
#' @param p2.upper upperbound for pseudotime2 (eg 1)
#' @param jmethod method for optimizing: use one that allows constraints like "L-BFGS-B"
#' @param jhessian return hessian (TRUE or FALSE). TRUE allows one to estimate the standard errors of the pseudotime estimates afterwards.
#' @return optim output for best pseudotime1 and pseudotime2 that fits the double-incubated cell
#' @examples
#' jfits.lst <- parallel::mclapply(jcells, function(jcell){
#'   cell.count.raw.merged <- mat.dbl[, jcell]
#'   jfit.fixedw <- FitMixingWeightPtimeLinear.fixedw(cell.count.raw.merged = cell.count.raw.merged, jfits.lowess1 = fits.lowess1.linear, jfits.lowess2 = fits.lowess2.linear, w.fixed = w.fixed, p1.init = 0.5, p1.lower = 0.01, p1.upper = 0.99,  p2.init = 0.5, p2.lower = 0.01, p2.upper = 0.99, jmethod = "L-BFGS-B", jhessian = TRUE)
#' }, mc.cores = ncores)
#' 
#' @export
FitMixingWeightPtimeLinear.fixedw <- function(cell.count.raw.merged, jfits.lowess1, jfits.lowess2, w.fixed, p1.init, p1.lower, p1.upper, p2.init, p2.lower, p2.upper, jmethod = "L-BFGS-B", jhessian = FALSE){
  # https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/hessian.nb.html
  # method by Brent is recommended?
  # jmethod <- ""L-BFGS-B""
  # w, ptime, cell.count.raw.merged, fits.lowess1, fits.lowess2
  optim(par = c(p1.init, p2.init),
        fn = GetLLMergedPtimeLowessLinear.fixedw,
        cell.count.raw.merged = cell.count.raw.merged,
        fits.lowess1 = jfits.lowess1,
        fits.lowess2 = jfits.lowess2,
        method = jmethod, hessian=jhessian, lower = c(p1.lower, p2.lower), upper = c(p1.upper, p2.upper), control = list(ndeps = c(0.02, 0.02)), w = w.fixed)
}
# list(fnscale = -1, trace=3)

