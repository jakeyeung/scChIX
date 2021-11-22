# Jake Yeung
# Date of Creation: 2020-04-03
# File: ~/projects/scchic-functions/scripts/processing_scripts/filter_cells_by_totalcuts_TAfrac_intrachromovar.R
# Stolen from filter_good_cells_good_bins.R, include intrachromovar. Remove blacklist filtering (should be done at bam level)
# 2020
# DESCRIPTION
# 
#     Filter cells by totalcuts_TAfrac_intrachromovar, optionaly remove EmptyWells
# 
# FOR HELP
# 
#     Rscript filter_cells_by_totalcuts_TAfrac_intrachromovar.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2019-12-24
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-infilerz', metavar='INFILESRZ', nargs = "+",
                                            help='infile containing TA frac')
parser$add_argument('-infilecounts', metavar='INFILESCOUNTMAT', nargs = "+",
                                            help='infile containing mat')
parser$add_argument('-names', metavar="space delim strings", nargs = "+",
                                            help='Label for each infile')
parser$add_argument('-outdir', metavar='OUTDIR',
                                            help='outdir')
parser$add_argument('-countcutoffmin', metavar='INTEGER', type = 'integer', nargs = "+", 
                                            help='Minimum counts')
parser$add_argument('-countcutoffmax', metavar='INTEGER', type = 'integer', nargs = "+", 
                                            help='Maximum counts')
parser$add_argument('-varcutoffmin', metavar='PosNumber', type = 'double', nargs = "+", 
                                            help='Minimum intrachromosomal variance')
parser$add_argument('-varcutoffmax', metavar='PosNumber', type = 'double', nargs = "+", 
                                            help='Maximum intrachromosomal variance')
parser$add_argument('-jmergesize', metavar='Nbins', type = 'integer', default=1000, 
                                            help='How many bins to merge when calculating global variance')
parser$add_argument('-TAcutoff', metavar='TAcutoff', type = 'double',
                                            help='Minimum TA fraction')
parser$add_argument('-chromoskeep', metavar='ChromosomeNames', type = 'character', nargs = "+", required=TRUE,
                                            help='List of chromosomes to keep')
parser$add_argument('--overwrite', action="store_true", default=FALSE, help="Force overwrite")
parser$add_argument('--doNotWriteTables', action="store_true", default=FALSE, help="Turn off writing rds outputs. Outputs pdf only")
parser$add_argument('--keepEmptyWells', action="store_true", default=FALSE, help="Do not removing the corner 8 wells we call EmptyWells")
parser$add_argument("--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")

                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

print(args)

print(args$names)
print(args$countcutoffmin)

# check files
assertthat::assert_that(length(args$infilerz) == length(args$infilecounts))
assertthat::assert_that(length(args$names) == length(args$infilecounts))
assertthat::assert_that(length(args$countcutoffmin) == 1 | length(args$countcutoffmin) == length(args$names), msg=paste("args$countcutoffmin length incompatible with args$names"))

# handles both length 1 or multillength 
cutoffmin.counts.hash <- hash::hash(args$names, args$countcutoffmin)
cutoffmax.counts.hash <- hash::hash(args$names, args$countcutoffmax)

varmin.counts.hash <- hash::hash(args$names, args$varcutoffmin)
varmax.counts.hash <- hash::hash(args$names, args$varcutoffmax)

print("making count hash:")
print(cutoffmin.counts.hash)

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

print("Starting time")
jstart <- Sys.time()

library(hash)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)

library(gtools)

# library(GenomicRanges)



# Constants set to arguments  ---------------------------------------------

outdir <- args$outdir
dir.create(outdir)

jmergesize <- args$jmergesize
cutoff.counts <- args$countcutoffmin
cutoff.TA <- args$TAcutoff
overwrite <- args$overwrite

infs.rz <- as.list(args$infilerz)
infs.mat <- as.list(args$infilecounts)

lapply(infs.rz, function(x) assertthat::assert_that(file.exists(x), msg=paste("infrz not found", x)))
lapply(infs.mat, function(x) assertthat::assert_that(file.exists(x), msg=paste("infmat not found", x)))

# jnames <- c("H3K4me1", "H3K27me3", "H3K4me1_H3K27me3")
jnames <- args$names
names(jnames) <- jnames

names(infs.rz) <- jnames
names(infs.mat) <- jnames

outdir <- args$outdir

outpaths <- lapply(jnames, function(jname){
  jout <- file.path(outdir, paste0("count_mat.",jname, ".countcutoffmin_", paste(cutoff.counts, collapse="-"), ".TAcutoff_", cutoff.TA, ".rds"))
  if (!overwrite){
    assertthat::assert_that(!file.exists(jout), msg = paste("Outfile exists, not overwriting for safety:", jout))
  }
  return(jout)
})

outbases <- lapply(jnames, function(jname){
  jout <- file.path(outdir, paste0("count_mat.",jname, ".countcutoffmin_", paste(cutoff.counts, collapse="-"), ".TAcutoff_", cutoff.TA))  # add extetnsion later 
  return(jout)
})

pdfout <- file.path(outdir, paste0("qc_plots.", paste(jnames, collapse = "-"), ".pdf"))
if (!overwrite){
  assertthat::assert_that(!file.exists(pdfout), msg = paste("Pdfout exists, not overwriting for safety:", pdfout))
}

names(infs.mat) <- jnames
names(infs.rz) <- jnames

print(infs.rz)

if (!any(sapply(infs.rz, function(inf) endsWith(inf, ".rds")))){
  dat.rz <- ReadLH.SummarizeTA(infs.rz, remove.nones = FALSE, na.to.zero = TRUE, bind.rows = FALSE)
} else {
  dat.rz <- lapply(infs.rz, function(inf) readRDS(inf))
}
# if (!any(endsWith(unlist(infs.rz), ".rds"))){
# } else {
# }

# add jname to dat.rz
dat.rz.filt <- lapply(jnames, function(jmark){
  jtmp <- dat.rz[[jmark]] %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(experi = ClipLast(samp, jsep = "_"),
         cellindx = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""))

# Load matrix -----------
# filter good cells
mats <- lapply(infs.mat, function(inf){
  if (endsWith(inf, ".csv") | endsWith(inf, ".txt")){
    print("Reading sliding win format")
    mat <- ReadMatSlideWinFormat(inf)
  } else if (endsWith(inf, ".rds")){
    print("Reading rds directly")
    mat <- readRDS(inf)
  } else {
    stop(inf, "must end with .csv or .rds")
  }
  # filter chromosomes
  jchromos.vec <- sapply(rownames(mat), function(x) strsplit(x, ":")[[1]][[1]])
  jchromos.filt.i <- which(jchromos.vec %in% args$chromoskeep)
  assertthat::assert_that(length(jchromos.filt.i) > 0)
  mat.filt <- mat[jchromos.filt.i, ]
  jchromos.check <- sort(unique(sapply(rownames(mat.filt), function(x) strsplit(x, ":")[[1]][[1]])))
  print("Chromos check")
  print(jchromos.check)

  # sort by rows in a reasonable way
  rorder <- gtools::mixedorder(rownames(mat.filt))
  mat.filt <- mat.filt[rorder, ]
  return(mat.filt)
})



# Get good cells  ---------------------------------------------------------

if (!args$keepEmptyWells){
  empty.wells <- GetEmptyWells()
} else {
  empty.wells <- c()  # empty
}

dat.rz.filt <- dat.rz.filt %>%
  rowwise() %>%
  mutate(empty.well = cellindx %in% empty.wells, 
         good.cell = total.count > cutoffmin.counts.hash[[mark]] & total.count < cutoffmax.counts.hash[[mark]] & TA.frac > cutoff.TA & !empty.well)


m.density.bymark <- ggplot(dat.rz.filt, aes(x = total.count, fill = mark)) + geom_density(alpha = 0.3)  + 
  scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_viridis_d()

m.scatter.bymark.col <- ggplot(dat.rz.filt, aes(x = total.count, y = TA.frac, color = good.cell, size = empty.well, shape = empty.well)) + geom_point()  + 
  scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1)

# for each mark, plot each plate 
m.scatter.lst <- lapply(jnames, function(jname){
  m <- ggplot(subset(dat.rz.filt, mark == jname), aes(x = total.count, y = TA.frac, color = good.cell, size = empty.well, shape = empty.well)) + geom_point()  +
    scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    facet_wrap(~mark, ncol = 1) + ggtitle(jname) + 
    geom_vline(xintercept = c(cutoffmin.counts.hash[[jname]], cutoffmax.counts.hash[[jname]])) + 
    geom_hline(yintercept = cutoff.TA)
  return(m)
})

m.density.lst <- lapply(jnames, function(jname){
  ggplot(subset(dat.rz.filt, mark == jname), aes(x = total.count, fill = mark)) + geom_density(alpha = 0.3)  +
    scale_x_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_viridis_d() + ggtitle(jname) + 
    geom_vline(xintercept = c(cutoffmin.counts.hash[[jname]], cutoffmax.counts.hash[[jname]]))
})


# how many good cells?
dat.rz.filt %>%
  group_by(experi) %>%
  filter(good.cell) %>%
  summarise(ncells = length(total.count))

dat.rz.filt %>%
  group_by(experi) %>%
  summarise(ncells = length(total.count))

cells.keep <- subset(dat.rz.filt, good.cell)$samp
print(paste("Number of cells before:", nrow(dat.rz.filt)))
print(paste("Number of cells after:", length(cells.keep)))




# Load mats ---------------------------------------------------------------

# filter good cells
mats <- lapply(mats, function(mat){
  cols.i <- colnames(mat) %in% cells.keep
  mat.filt <- mat[, cols.i]
  return(mat.filt)
})



# Load variance and calculate it  -----------------------------------------

dat.vars.raw.lst <- lapply(jnames, function(jname){
  dat.var <- CalculateVarRaw(mats[[jname]], merge.size = jmergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"))
  dat.var$mark <- jname
  dat.var <- dat.var %>%
    rowwise() %>%
    mutate(good.cells.var = ncuts.var > varmin.counts.hash[[jname]] & ncuts.var < varmax.counts.hash[[jname]])
  return(dat.var)
})


# Plot the var  -----------------------------------------------------------

m.var.platesmerged.lst <- lapply(jnames, function(jname){
  dat.vars.raw <- dat.vars.raw.lst[[jname]]
  m.var <- ggplot(dat.vars.raw, aes(x = ncuts, y = ncuts.var, color = good.cells.var)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_log10() + scale_y_log10()  +  ggtitle(jname) + 
    geom_hline(yintercept = c(varmin.counts.hash[[jname]],  varmax.counts.hash[[jname]]))
 return(m.var) 
})

m.var.lst <- lapply(jnames, function(jname){
  dat.vars.raw <- dat.vars.raw.lst[[jname]]
  m.var <- ggplot(dat.vars.raw, aes(x = ncuts, y = ncuts.var, color = good.cells.var)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_log10() + scale_y_log10()  +  ggtitle(jname) + 
    facet_wrap(~plate) + geom_hline(yintercept = c(varmin.counts.hash[[jname]],  varmax.counts.hash[[jname]]))
 return(m.var) 
})

# label badvarcells 
dat.vars.raw.long.filt <- dat.vars.raw.lst %>%
  bind_rows() %>%  
  filter(good.cells.var)


cells.keep.goodvar <- dat.vars.raw.long.filt$cell


# Filter by var -----------------------------------------------------------

print("Dimensions before filtering by var:")
print(lapply(mats, dim))
# filter good cells
mats <- lapply(mats, function(mat){
  cols.i <- colnames(mat) %in% cells.keep.goodvar
  mat.filt <- mat[, cols.i]
  return(mat.filt)
})
print("Dimensions after filtering by var:")
print(lapply(mats, dim))


# Plot pdf outputs --------------------------------------------------------

pdf(pdfout, width = 1440/72, height = 815/72, useDingbats = FALSE)
  print(m.scatter.bymark.col)
  print(m.density.bymark)
  print(m.scatter.lst)
  print(m.density.lst)
  print(m.var.lst)
  print(m.var.platesmerged.lst)
dev.off()

if (!args$doNotWriteTables){
  # save each mat separately 
  for (jname in jnames){
    outpath <- outpaths[[jname]]
    jtmp <- mats[[jname]]
    assertthat::assert_that(nrow(jtmp) > 0 & ncol(jtmp) > 0)
    saveRDS(jtmp, file = outpath)
    writeMM(jtmp, sparse = TRUE, file = paste0(outbases[[jname]], ".mm"))
    write.table(rownames(jtmp), file = paste0(outbases[[jname]], ".rownames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    write.table(colnames(jtmp), file = paste0(outbases[[jname]], ".colnames"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
  }
} else {
  print(paste("Do not write tables:", args$doNotWriteTables))
  print("Skipping writing, but plots still written")
}

