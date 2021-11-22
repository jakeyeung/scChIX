# Jake Yeung
# run_LDA_model2.R
# 2019-11-10
# DESCRIPTION
# 
#     Run LDA clean up params
# 
# FOR HELP
# 
#     Rscript run_LDA_model2.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-11-10
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

jstart <- Sys.time() 
suppressPackageStartupMessages(library("argparse"))
library(topicmodels)
library(dplyr)
library(ggplot2)
library(Matrix)
library(here)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument('inpath', metavar='INFILE',
                                            help='.RData where count mat is in count.dat$counts or .rds object to count mat')
parser$add_argument('outpath', metavar='OUTFILE',
                                            help='Output .Robj of LDA outputs')
parser$add_argument("-t", "--topics", metavar='Comma sep string', required=TRUE,
                                            help='CSV of topics to iterate')
parser$add_argument("-b", "--binarizemat", action="store_true", default=FALSE,
                        help="Binarize matrix")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
parser$add_argument("--RemoveDupRows", action="store_true", default=FALSE,
                        help="Remove duplicated rows, default FALSE")
parser$add_argument("--RemoveEmptyCells", action="store_true", default=FALSE,
                        help="Remove empty cols, default FALSE")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()


setwd(here())
print(paste("Work directory: ", getwd()))
# source("scripts/Rfunctions/ParseStrings.R")
# source("scripts/Rfunctions/Aux.R")

# args <- commandArgs(trailingOnly=TRUE)

print("input args:")
print(args)

inpath <- args$inpath
outpath <- args$outpath
topic.vec <- as.numeric(StrToVector(args$topics, delim = ","))
binarizemat <- args$binarizemat

print(paste("Will iterate through", length(topic.vec), "Ks"))
print(topic.vec)


# Load counts -------------------------------------------------------------

if (endsWith(inpath, ".rds")){
  count.mat <- readRDS(inpath)
} else {
  # assume .RData object with count.dat$counts
  load(inpath, v=T)
  count.mat <- count.dat$counts
}

# remove empty cols
if (args$RemoveEmptyCells){
  print("Removing empty cells...")
  print(dim(count.mat))
  cols.empty <- colSums(count.mat) == 0
  count.mat <- count.mat[, !cols.empty]
  print(dim(count.mat))
}

print(dim(count.mat))

if (args$RemoveDupRows){
  print("Removing dup rows")
  rnames.all <- rownames(count.mat)
  rnames.keep <- !duplicated(rnames.all)
  dim(count.mat)
  count.mat <- count.mat[rnames.keep, ]
  print("Mat after removing dup rows")
  dim(count.mat)

  print("Removing dup cols")
  cnames.all <- colnames(count.mat)
  cnames.keep <- !duplicated(cnames.all)
  dim(count.mat)
  count.mat <- count.mat[, cnames.keep]
  print("Mat after removing dup cols")
  dim(count.mat)
}


# Run LDA on count matrix -------------------------------------------------

print("Running LDA")

# binarize matrix
if (binarizemat){
  count.mat.orig <- count.mat
  count.mat <- BinarizeMatrix(count.mat)
  print(paste('Max count after binarizing', max(count.mat)))
} else {
  count.mat.orig <- NA
}

if (length(topic.vec) > 1){
  print("Running multicore LDA for topics:")
  print(topic.vec)
  out.lda <- parallel::mclapply(topic.vec, function(nc){
          LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))            
  }, mc.cores = length(topic.vec))
} else {
  print("Running single LDA for topics:")
  print(topic.vec)
  out.lda <- LDA(x = t(count.mat), k = topic.vec, method = "Gibbs", control=list(seed=0))
}

# save output
print("Saving LDA")
save(out.lda, count.mat, count.mat.orig, file = outpath)
print("Time elapsed after LDA")
print(Sys.time() - jstart)



print("Time elapsed after tuning")
print(Sys.time() - jstart)

