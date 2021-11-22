# Jake Yeung
# Date of Creation: 2020-02-05
# File: ~/projects/dblchic/scripts/macbook_analysiis/pretty_analysis_EtOH/setup_objs_for_unmixing_topics_general.R
# General with arguments


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)


# Contants ----------------------------------------------------------------

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument('-mark1', metavar='STR',
                    help='First mark (e.g. active mark)')
parser$add_argument('-mark2', metavar='STR',
                    help='Second mark (e.g., repressive mark)')
parser$add_argument('-inf1', metavar='PATH',
                    help='mark1 Path to UMAP annotated data containing LDA output (out.lda) act and repress as well as celltype annotations for each cluster (dat.merge, with cluster colname)')
parser$add_argument('-inf2', metavar='PATH',
                    help='mark2 Path to UMAP annotated data containing LDA output (out.lda) for act and repress as well as celltype annotations for each cluster (dat.merge, with cluster colname)')
parser$add_argument('-infdbl', metavar='PATH',
                    help='dblmark Path to UMAP annotated data containing LDA output (out.lda, count.mat) act and repress as well as celltype annotations for each cluster (dat.merge, with cluster colname, for removing NAs). count.mat used to set up raw counts')
parser$add_argument('-outprefix', metavar='OUTFILE',
                    help='Prefix to write .RData and .pdf')
parser$add_argument("--RemoveNA", action="store_true", default=FALSE,
                    help="Any cluster that is assigned NA will be excluded as a cluster and further analyses. If FALSE, use NAs as separate cluster and keep raw dbl marks marked as NA for unmixing")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()


remove.na <- args$RemoveNA
jmark <- args$mark1
jmark2 <- args$mark2
jmarks <- c(jmark, jmark2); names(jmarks) <- jmarks

jmarks.str <- paste(jmark, jmark2, sep = "-")

jmark.active <- jmarks[[1]]
jmark.repress <- jmarks[[2]]

# init output file 
outfprefix <- args$outprefix
# add RemoveNA TRUE or FALSE in outfprefix
outfprefix <- paste0(outfprefix, ".removeNA_", remove.na)
outfobjs <- paste0(outfprefix, ".RData")
outpdf <- paste0(outfprefix, ".pdf")
assertthat::assert_that(!file.exists(outfobjs))


inf1 <- args$inf1
assertthat::assert_that(file.exists(inf1))
inf2 <- args$inf2
assertthat::assert_that(file.exists(inf2))
inf.dbl <- args$infdbl
assertthat::assert_that(file.exists(inf.dbl))

# UMAP settings -----------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# colors for UMAP 
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# Load LDA outputs --------------------------------------------------------


load(inf1, v=T)
tm.result.act <- posterior(out.lda)
dat.merge.act <- dat.merge

load(inf2, v=T)
tm.result.repress <- posterior(out.lda)
dat.merge.repress <- dat.merge

tm.result.lst <- list(tm.result.act, tm.result.repress)
names(tm.result.lst) <- jmarks

dat.merge.lst <- list(dat.merge.act, dat.merge.repress)
names(dat.merge.lst) <- jmarks

dat.louv <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.louv <- dat.merge.lst[[jmark]]
  if (remove.na){
    dat.louv <- subset(dat.louv, !is.na(cluster))
  }
  dat.louv$mark <- jmark
  dat.louv$cluster <- gsub(pattern = "_", replacement = "", dat.louv$cluster)  # alow only one underscore  # cal lit CLUSTER
  return(dat.louv)
})


mlst <- lapply(jmarks, function(jmark){
  x <- dat.louv[[jmark]]
  m <- ggplot(x, aes(x = umap1, y = umap2, color = cluster)) + geom_point(alpha = 1, size = 3) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark)
})

# write to output file
pdf(outpdf, useDingbats = FALSE)
  print(mlst)
dev.off()

# Get raw  ----------------------------------------------------------------

load(inf.dbl, v=T)
# filter count.mat for non NA cells
if (remove.na){
  cells.keep <- unique(subset(dat.merge, !is.na(cluster))$cell)  # unique prevents cases where dat.merge contains duplicated cells!
  print("Dim raw counts before filtering NAs...")
  print(dim(count.mat))
  count.mat <- count.mat[, cells.keep]
  print("Dim raw counts after filtering NAs...")
  print(dim(count.mat))
} else {
  cells.keep <- dat.merge$cell
  print("Dim raw counts, no filtering of NAs")
  print(dim(count.mat))
}

rnames.lst <- list(colnames(tm.result.lst[[jmark.active]]$terms), colnames(tm.result.lst[[jmark.repress]]$terms), rownames(count.mat))
rnames.common <- Reduce(f = intersect, x = rnames.lst)

print(length(rnames.common))
# # get common rows?
# rnames.common <- intersect(rownames(count.mat), colnames(dat.impute.active))
count.dat <- list()
count.dat$counts <- count.mat[rnames.common, ]

# Take average cluster as my "topics" -------------------------------------

dat.imputes <- lapply(jmarks, function(jmark){
  dat.impute <- tm.result.lst[[jmark]]$topics %*% tm.result.lst[[jmark]]$terms
  dat.impute <- dat.impute[, rnames.common]
  print(dim(dat.impute))
  return(dat.impute)
})

# average across cells
cell.avgs <- lapply(jmarks, function(jmark){
  jsplit <- split(dat.louv[[jmark]], dat.louv[[jmark]]$cluster)
  return(jsplit)
})


# Get cell-averaged proportions for each louvain  -------------------------

dat.impute.repress.lst <- lapply(cell.avgs[[jmark.repress]], function(jsub.louvain){
  jcells <- as.character(jsub.louvain$cell)
  jrow <- colMeans(dat.imputes[[jmark.repress]][jcells, ])
  return(jrow)
})

# make dat.active.mat: columns are genomic regions, rows are louvains?
dat.impute.active <- lapply(cell.avgs[[jmark.active]], function(jsub.louvain){
  jcells <- as.character(jsub.louvain$cell)
  jrow <- colMeans(dat.imputes[[jmark.active]][jcells, ])
}) %>%
  bind_rows() %>%
  as.data.frame() %>%
  t()

colnames(dat.impute.active) <- colnames(dat.imputes[[jmark.active]])

# check names
assertthat::assert_that(all(colnames(dat.impute.active) == names(dat.impute.repress.lst[[1]])))
assertthat::assert_that(all(rownames(count.dat$counts) == names(dat.impute.repress.lst[[1]])))

save(count.dat, dat.impute.repress.lst, dat.impute.active, dat.louv, file = outfobjs)

