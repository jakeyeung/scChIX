# Jake Yeung
# Date of Creation: 2019-08-05
# Double staining functions


#' Given a mixing fraction w, two probability vectors, split counts from double-incubated reads into their respective histone modification.
#' 
#' @param x.raw vector of raw counts from a double-incubated cell.
#' @param mixweight mixing fraction of the two histone modifications, relative to a reference histone modification (usually p.active)
#' @param p.active probability vector across genomic bins for reference histone modification. ie probability a sampled read from reference histone modification falls on specific genomic locus.
#' @param p.repress probability vector across genomic bins for non-reference histone modification. ie probability a sampled read from non-reference histone modification falls on specific genomic locus.
#' @param random.seed random seed for reproducibility
#' @return list: x.raw.active: a vector of counts assigned to reference histone mark. x.raw.repress: a vector of counts assigned to non-ref histone mark, p.cell.active.weights: parameter used to specify the binomial distirbution for splitting reads.
#' @examples
#' x.raw.unmixed <- lapply(all.cells, function(jcell){
#'   # print(jcell)
#'   return(UnmixRawCounts(x.raw = all.x.raw[[jcell]], mixweight = all.mixweights[[jcell]], p.active = all.p.active[[jcell]], p.repress = all.p.repress[[jcell]], random.seed = 0))
#' })
#' )
#' @export
UnmixRawCounts <- function(x.raw, mixweight, p.active, p.repress, random.seed = 0){
  # unmix cells into active and repress using binomial
  set.seed(random.seed)
  p.cell.active.weights <- mixweight * p.active /  ( mixweight * p.active + (1 - mixweight) * p.repress)
  assertthat::assert_that(length(p.cell.active.weights) == length(x.raw))
  # p.cell.repress.weights <- (1 - mixweight) * p.repress /  (mixweight * p.active + (1 - mixweight) * p.repress)
  x.raw.active <- rbinom(n = length(x.raw), size = x.raw, prob = p.cell.active.weights)
  x.raw.repress <- x.raw - x.raw.active  # deterministic
  return(list(x.raw.active = x.raw.active, x.raw.repress = x.raw.repress, p.cell.active.weights = p.cell.active.weights))
}


GetDoubleName <- function(jname.test, cell.lambda.hash, prefix = "Pseudo:"){
  # used for simulation
  return(paste0(prefix, paste(lapply(strsplit(jname.test, ";")[[1]], function(x) cell.lambda.hash[[x]]), collapse = ";")))
}

#' Calculate log-likelihood for all pairs of clusters
#' 
#' @param w mixing fraction of the two histone modifications, relative to a reference histone modification (usually dat.impute.active)
#' @param cell.count.raw.merged vector of raw cell counts from a double-incubated cell. 
#' @param dat.impute.repress.lst Non-reference histone modification: list object containing cluster-specific probability weights across genomic bins. Each element in list is a probability vector across genomic bins for a cluster.
#' @param dat.impute.active Reference histone modification: matrix object where rows are cluster names and columns are genomic bins. Each row is a probability vector across genomic bins. 
#' @param return.mat if TRUE, returns full log-likelihood matrix for all possible pairs. If False, returns the negative of the maximum likelihood value.
#' @return if return.mat=TRUE, returns full log-likelihood matrix for all possible pairs. If False, returns the negative of the maximum likelihood value.
#' @examples
#' data(RawDblCountMatSubset)
#' cell.count.raw.merged.lst <- as.list(as.data.frame(as.matrix(count.mat.dbl.subset)))
#' act.repress.coord.lst <- lapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
#'   optim.out <- FitMixingWeight(cell.count.raw.merged = cell.count.raw.merged,
#'                                dat.impute.repress.lst = dat.impute.repress.lst,
#'                                dat.impute.active = dat.impute.active, w.init = 0.5, w.lower = wlower, w.upper = wupper, jmethod = "Brent")
#'   ll.mat <- GetLLMerged(optim.out$par, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = TRUE)
#'   return(list(ll.mat = ll.mat, w = optim.out$par))
#' })
#' @export
GetLLMerged <- function(w, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = FALSE){
  w.vec <- c(w, 1-w)  # active proportion, repress proportion
  if (any(w.vec < 0)){
    print(paste("Error message: w is negative:"))
    print(w.vec)
    w.vec[which.min(w.vec)] <- 0
    w.vec[which.max(w.vec)] <- 1
    print("Fixed w:")
    print(w.vec)
  }
  dat.imputed.merged <- lapply(dat.impute.repress.lst, function(repress.row){
    t(outmat <- apply(dat.impute.active, 1, function(act.row){
      matrixStats::colWeightedMeans(rbind(act.row, repress.row), w = w.vec)  # standard weighted.mean is SUPER slow
    }))
  })
  # Calculate all the likelihoods  ------------------------------------------
  # output will be N active cell by N repress cell
  repress.names <- names(dat.imputed.merged)
  names(repress.names) <- repress.names

  ll.lst <- lapply(repress.names, function(repress.name){
    repress.mat <- dat.imputed.merged[[repress.name]]
    ll <- as.data.frame(apply(repress.mat, 1, function(merged.prob.vec){
      return(dmultinom(x = cell.count.raw.merged, prob = merged.prob.vec, log = TRUE))
    }))
    rownames(ll) <- rownames(repress.mat)
    colnames(ll) <- repress.name
    return(ll)
    # return(repress.mat)
  })
  ll.merged <- do.call(cbind, ll.lst)
  if (return.mat){
    return(ll.merged)
  } else {
    # return negative of maximum likelihood
    # negative because optimizers find minimums
    return(-1 * max(ll.merged))
  }
}


#' Infer best pair of clusters from a double-incubated cell that maximizes the multinomial log-likelihood
#' 
#' @param cell.count.raw.merged vector of raw cell counts from a double-incubated cell. 
#' @param dat.impute.repress.lst Non-reference histone modification: list object containing cluster-specific probability weights across genomic bins. Each element in list is a probability vector across genomic bins for a cluster.
#' @param dat.impute.active Reference histone modification: matrix object where rows are cluster names and columns are genomic bins. Each row is a probability vector across genomic bins. 
#' @param w.init initial guess for the mixing fraction between the two histone modifications. Mixing fraction is relative to the reference histone modification.
#' @param w.lower lower bound for mixing weights, usually 0.
#' @param w.upper lower bound for mixing weights, usually 1.
#' @param jmethod method for optimizing w. Use "Brent" because it is a 1D problem given two clusters.
#' @return optim output for best pair of clusters (one from an element in dat.impute.repress.lst, other from a row from dat.impute.active).
#' @examples
#' data(RawDblCountMatSubset)
#' cell.count.raw.merged.lst <- as.list(as.data.frame(as.matrix(count.mat.dbl.subset)))
#' act.repress.coord.lst <- lapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
#'   optim.out <- FitMixingWeight(cell.count.raw.merged = cell.count.raw.merged,
#'                                dat.impute.repress.lst = dat.impute.repress.lst,
#'                                dat.impute.active = dat.impute.active, w.init = 0.5, w.lower = wlower, w.upper = wupper, jmethod = "Brent")
#'   ll.mat <- GetLLMerged(optim.out$par, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = TRUE)
#'   return(list(ll.mat = ll.mat, w = optim.out$par))
#' })
#' @export
FitMixingWeight <- function(cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, w.init, w.lower, w.upper, jmethod = "Brent"){
  # https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/hessian.nb.html
  # method by Brent is recommended?
  # jmethod <- ""L-BFGS-B""
  optim(par = w.init,
        fn = GetLLMerged,
        cell.count.raw.merged = cell.count.raw.merged,
        dat.impute.repress.lst = dat.impute.repress.lst,
        dat.impute.active = dat.impute.active,
        method = jmethod, hessian=TRUE, lower = w.lower, upper = w.upper)
}
