# Jake Yeung
# load_objs_run_fits.R
# 2019-08-20
# DESCRIPTION
# 
#     Run objects for double staining and run fits
# 
# FOR HELP
# 
#     Rscript load_objs_run_fits.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-08-20
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(scchicFuncs)
library(umap)
library(topicmodels)
library(irlba)
library(hash)
library(igraph)

library(Matrix)
library(parallel)
library(scchicFuncs)

print(getwd())

library(here)
setwd(here())
print("Current wd")
print(getwd())
# source("scripts/Rfunctions/DoubleStaining.R")

parser$add_argument('infile', metavar='INFILE',
                                            help='Input .RData containing count.dat, dat.impute.repress.lst, dat.impute.active, dat.louv')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output from fits')
parser$add_argument("-n", "--ncores", type='integer', default=4,
                        help="Number of cores to run (default 4)")
parser$add_argument("-m", "--method", type='character', default="Brent",
                        help="Method for optimization")
parser$add_argument("-wlower", type='double', default=0,
                        help="Minimum w, mixing fraction for first mark")
parser$add_argument("-wupper", type='double', default=1,
                        help="Maximum w, mixing fraction for first mark")
parser$add_argument("-pseudocount", type='integer', default=0,
                        help="Pseudocount")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

inf.objs <- args$infile
outf <- args$outfile
assertthat::assert_that(file.exists(inf.objs))

load(inf.objs, v=T)

# get list of cell.counet.raw.merged
cell.count.raw.merged.lst <- as.list(as.data.frame(as.matrix(count.dat$counts)))
system.time(
  act.repress.coord.lst <- mclapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
    if (args$pseudocount == 0){
        # print("No pseudocount added")
    } else {
        print(paste("Adding pseudocount", args$pseudocount))
        cell.count.raw.merged <- cell.count.raw.merged + args$pseudocount
    }
    optim.out <- FitMixingWeight(cell.count.raw.merged = cell.count.raw.merged, 
                                 dat.impute.repress.lst = dat.impute2.lst, 
                                 dat.impute.active = dat.impute1, w.init = 0.5, w.lower = args$wlower, w.upper = args$wupper, jmethod = args$method)
    ll.mat <- GetLLMerged(optim.out$par, cell.count.raw.merged, dat.impute2.lst, dat.impute1, return.mat = TRUE)
    return(list(ll.mat = ll.mat, w = optim.out$par))
  }, mc.cores = args$ncores)
)

save(act.repress.coord.lst, file = outf)
