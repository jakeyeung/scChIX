# Jake Yeung
# project_new_samples_on_LDA.R
# 2019-06-18
# DESCRIPTION
# 
#     Project new samples (in form of a sparse matrix) onto existing LDA
# 
# FOR HELP
# 
#     Rscript project_new_samples_on_LDA.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-06-18
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(topicmodels)
library(JFuncs)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE RDATA',
                                            help='Input .RData containing LDA object named out.objs$out.lda or out.lda. Or .rds of LDA_Gibbs class')
parser$add_argument('inmat', metavar='INMAT RDS',
                                            help='Input .rds containing sparse matrix (bins in rows, samples in columns). Or .RData of count.dat list object with count.dat$counts slot. dgCMatirx class')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Output RData containing posterior of predicted counts')
parser$add_argument("-b", "--binarizemat", action="store_true", default=FALSE,
                        help="Binarize matrix")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE,
                        help="Print extra output [default]")
parser$add_argument("--RemoveEmptyCells", action="store_true", default=FALSE,
                        help="Remove cells with zero reads. Means mixings were either 0 or 1")
parser$add_argument("--MatchRowsPadZeros", action="store_true", default=FALSE,
                        help="Match rows and pad zeros for rows that are missing")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

if (endsWith(args$infile, ".rds")){
    out.lda <- readRDS(args$infile)
    out.objs <- list()
    out.objs$out.lda <- out.lda
} else {
    objnames <- load(args$infile, v=T)  # out.objs$out.lda or out.lda
    if (any(objnames == "out.objs")){
        # dont do any rename 
        print("Object loaded is out.objs, doing nothing")
    } else if (any(objnames == "out.lda")){
        # rename
        print("Object loaded is out.lda, renaming to out.objs with slot out.lda")
        out.objs <- list()
        out.objs$out.lda <- out.lda
    }
}

assertthat::assert_that(class(out.objs$out.lda) == "LDA_Gibbs")

if (endsWith(args$inmat, ".rds")){
    print("Reading .rds file into count.dat") 
    count.dat <- list()
    count.dat$counts <- readRDS(args$inmat)
} else {
    print("Assuming .RData list object named count.dat with counts slot") 
    load(args$inmat, v=T)  # count.dat$counts
}

if (args$RemoveEmptyCells){
  # check size of cells, remove ones that are zero??
  print("Number of cells before filtering:")
  print(ncol(count.dat$counts))
  cells.keep <- which(colSums(count.dat$counts) > 0)
  print("Number of cells after filtering:")
  count.dat$counts <- count.dat$counts[, cells.keep]
  print(ncol(count.dat$counts))
}
assertthat::assert_that(class(count.dat$counts) == "dgCMatrix" | class(count.dat$counts) == "matrix")

if (args$binarizemat){
  print(paste('Max count before binarizing', max(count.dat$counts)))
  count.dat$counts <- BinarizeMatrix(count.dat$counts)
  print(paste('Max count after binarizing', max(count.dat$counts)))
}

if (args$MatchRowsPadZeros){

  print("Nterms mat dim:")
  print(length(out.objs$out.lda@terms))

  print("Input mat dim before:")
  print(dim(count.dat$counts))

  rnames.orig <- out.objs$out.lda@terms

  rows.filt1 <- rownames(count.dat$counts) %in% rnames.orig
  mat.filt.forproj1 <- count.dat$counts[rows.filt1, ]

  # # Match rows 
  # rnames.common <- intersect(rnames.orig, rownames(count.dat$counts))
  # assertthat::assert_that(length(rnames.common) > 0)
  # mat.filt.forproj1 <- count.dat$counts[rnames.common, ]
  
  # pad zeros
  rnames.add <- rnames.orig[!rnames.orig %in% rownames(count.dat$counts)]
  mat.add <- matrix(data = 0, nrow = length(rnames.add), ncol = ncol(mat.filt.forproj1), dimnames = list(rnames.add, colnames(mat.filt.forproj1)))
  mat.filt.forproj <- Matrix(rbind(mat.filt.forproj1, mat.add), sparse = TRUE)
  # reearrange rows
  mat.filt.forproj <- mat.filt.forproj[rnames.orig, ]
  count.dat$counts <- mat.filt.forproj
  print("Input mat dim after:")
  print(dim(count.dat$counts))

  print("Range of values input mat:")
  print(range(count.dat$counts))
}

count.mat.proj <- count.dat$counts


print("Projecting new samples onto trained LDA")
system.time(
  out.lda.predict <- posterior(out.objs$out.lda, t(as.matrix(count.dat$counts)))
)

save(out.lda.predict, count.mat.proj, out.objs, file = args$outfile)

