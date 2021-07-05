# Jake Yeung
# rbind_RZ_files.R
# 2021-06-28
# DESCRIPTION
# 
#     Rbind RZ files across replicates
# 
# FOR HELP
# 
#     Rscript rbind_RZ_files.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2021-06-28
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))
library(scchicFuncs)
library(data.table)
library(dplyr)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infile', metavar='INFILE', nargs = "+", 
                                            help='Space delim path to RZ files')
parser$add_argument('-outfile', metavar='OUTFILE',
                                            help='Output rds')
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

dat.rz <- ReadLH.SummarizeTA(args$infile, remove.nones = FALSE, na.to.zero = TRUE, bind.rows = TRUE)
print(dim(dat.rz))
saveRDS(dat.rz, file = args$outfile)
