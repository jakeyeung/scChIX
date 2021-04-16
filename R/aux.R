# Jake Yeung
# Date of Creation: 2021-04-16
# File: ~/projects/scChIX/R/aux.R
#


ClipLast <- function(x, jsep = "-", jsep.out = NULL){
  # B6-13W1-BM-H3K4me3-1_269 -> B6-13W1-BM-H3K4me3
  if (is.null(jsep.out)){
    jsep.out <- jsep
  }
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit) - 1
  return(paste(jsplit[1:N], collapse = jsep.out))
}



SoftMax <- function(x, return.log = TRUE, logfn = log){
  # numericallys table softmax by subtracting the maximum first
  # https://stackoverflow.com/questions/42599498/numercially-stable-softmax
  # x value are in log if .log = TRUE
  # numer <- log(exp(x - max(x)))
  # denom <- log(sum(exp(x - max(x))))
  numer <- logfn(exp(x - max(x)))
  denom <- logfn(sum(exp(x - max(x))))
  plog <- numer - denom
  if (return.log){
    return(plog)
  } else {
    return(exp(plog))
  }
}

