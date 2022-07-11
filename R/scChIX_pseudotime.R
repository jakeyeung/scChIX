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

PredictSignalGenomeWideLowessLinear <- function(lowess.fits, ptime, linearize = TRUE, normalize = TRUE){
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

