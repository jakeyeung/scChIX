# Jake Yeung
# Date of Creation: 2020-02-06
# File: ~/projects/dblchic/scripts/macbook_analysiis/pretty_analysis_EtOH/unmixing_downstream_split_reads_for_server.R
# Unmixing downstream general enough to be used as a script
# Handle downstream if we take union of rows then lda inputs don't match raw double incubated counts

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
library(JFuncs)

library(irlba)


# Parse args --------------------------------------------------------------

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument('-mark1', metavar='STR', required = TRUE,
                    help='First mark (e.g. active mark)')
parser$add_argument('-mark2', metavar='STR', required = TRUE,
                    help='Second mark (e.g., repressive mark)')
parser$add_argument('-inf_dbl_input', metavar='PATH to RDATA', required = TRUE,
                    help='Input .RData for double staining inference')
parser$add_argument('-inf_dbl_output', metavar='PATH to RDATA', required = TRUE,
                    help='Output .RData for double staining inference')
parser$add_argument('-inf_mark1', metavar='PATH to RDATA or ROBJ', required = TRUE,
                    help='mark1 outpath containig count.mat or RDS of CountMat')
parser$add_argument('-inf_mark2', metavar='PATH to RDATA or ROBJ', required = TRUE,
                    help='mark2 outpath containing count.mat or RDS of CountMat')
parser$add_argument('-inf_dblmark_lda', metavar='PATH to RDATA or ROBJ', required = TRUE,
                    help='Double mark LDA output path')
# parser$add_argument('-inf_mark1_countmat', metavar='PATH to RDATA or ROBJ', required = TRUE,
#                     help='mark1 CountMat output path')
# parser$add_argument('-inf_mark2_lda', metavar='PATH to RDATA or ROBJ', required = TRUE,
#                     help='mark2 LDA output path')
# parser$add_argument('-inf_dblmark_countmat', metavar='PATH to RDATA or ROBJ', required = TRUE,
#                     help='Double mark countmat output path')
parser$add_argument('-outprefix', metavar='OUTFILE', required = TRUE,
                    help='Prefix to write .RData and .pdf')
parser$add_argument('--RemoveBadMixings', default=FALSE, action="store_true",
                    help='Remove cells with mixings less than 0.01 or greater than 0.99')


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()


# Set constants -----------------------------------------------------------



jcutoff <- 0.95
logpcutoff <- log(jcutoff)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jtxtsize <- 4
# cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")


# Set variables from args -------------------------------------------------


jmark1 <- args$mark1
jmark2 <- args$mark2

outprefix <- args$outprefix

inf.in <- args$inf_dbl_input
assertthat::assert_that(file.exists(inf.in))
inf.out <- args$inf_dbl_output
assertthat::assert_that(file.exists(inf.out))
# dbl also needed later for splitting reads

inf.act <- args$inf_mark1
assertthat::assert_that(file.exists(inf.act))

inf.rep <- args$inf_mark2
assertthat::assert_that(file.exists(inf.act))
assertthat::assert_that(file.exists(inf.rep))

inf.dbl.lda <- args$inf_dblmark_lda
assertthat::assert_that(file.exists(inf.dbl.lda))

# inf.act.mat <- args$inf_mark1_countmat
# inf.rep.mat <- args$inf_mark2_countmat
# inf.dbl.mat <- args$inf_dblmark_countmat
# assertthat::assert_that(file.exists(inf.act.mat))
# assertthat::assert_that(file.exists(inf.rep.mat))
# assertthat::assert_that(file.exists(inf.dbl.mat))

# Set paths ---------------------------------------------------------------

jmarks <- c(jmark1, jmark2); names(jmarks) <- jmarks
jmarks.pair <- c(jmark1, jmark2)
jmark.dbl.dash <- paste(jmarks.pair, collapse = "-")
# jmark.active <- jmarks[[1]]
# jmark.repress <- jmarks[[2]]
# jmark.act <- jmarks.pair[[1]]
# jmark.repress <- jmarks.pair[[2]]
jmark.dbl <- paste(c(jmark1, jmark2), collapse = "-")
names(jmarks.pair) <- jmarks.pair

outpdf <- paste0(outprefix, "-plots.pdf")
outpdf.mats <- paste0(outprefix, "-plots_matsLL.pdf")
outpdf.pairs <- paste0(outprefix, "-plots_pairs.pdf")

# outputs for splitting reds 
outpdf.split <- paste0(outprefix, "-plots_WithUnmixing.pdf")
outf.unmixing <- paste0(outprefix, "-raw_unmixed.RData")
outf.unmixing.probs.mat <- paste0(outprefix, "-prob_mat.", jmark.dbl.dash, "_to_", jmark1, ".txt")

outf.mergedmat.act <- paste0(outprefix, "-merged_mat.", jmark1, ".rds")
outf.mergedmat.repress <- paste0(outprefix, "-merged_mat.", jmark2, ".rds")
outf.unmixedmat.act <- paste0(outprefix, "-unmixed_mat.", jmark1, ".rds")
outf.unmixedmat.repress <- paste0(outprefix, "-unmixed_mat.", jmark2, ".rds")

outtxt.badcells <- paste0(outprefix, "-badcells.txt")
outrdata.badcells <- paste0(outprefix, "-badcells.RData")


# Filenames ---------------------------------------------------------------




# Load double  ------------------------------------------------------------


load(inf.in, v=T)
count.mat.dbl <- count.dat$counts
load(inf.out, v=T)
rnames.keep <- rownames(count.mat.dbl)

# 2020-10-22
# rnames.keep has more rows than inf.act and inf.rep, need to handle this properly 

if (endsWith(inf.act, "rds")){
  count.mat.act.orig <- readRDS(inf.act)
} else {
  load(inf.act, v=T)
  count.mat.act.orig <- count.mat
}

rnames.act <- rownames(count.mat.act.orig)  # for later when we force count.mat.dbl to be same rownames as LDA object terms
rnames.add.act <- rnames.keep[which(!rnames.keep %in% rownames(count.mat.act.orig))]
mat.add.act <- matrix(0, nrow = length(rnames.add.act), ncol = ncol(count.mat.act.orig), dimnames = list(rnames.add.act, colnames(count.mat.act.orig)))

count.mat.act <- rbind(count.mat.act.orig, mat.add.act)
count.mat.act <- count.mat.act[rnames.keep, ]

if (endsWith(inf.rep, ".rds")){
  count.mat.repress.orig <- readRDS(inf.rep)
} else {
  load(inf.rep, v=T)
  count.mat.repress.orig <- count.mat
}

rnames.repress <- rownames(count.mat.repress.orig)  # to make count.mat.dbl same rownames as original count.mat.rep

rnames.add.repress <- rnames.keep[which(!rnames.keep %in% rownames(count.mat.repress.orig))]
mat.add.repress <- matrix(0, nrow = length(rnames.add.repress), ncol = ncol(count.mat.repress.orig), dimnames = list(rnames.add.repress, colnames(count.mat.repress.orig)))

count.mat.repress <- rbind(count.mat.repress.orig, mat.add.repress)
count.mat.repress <- count.mat.repress[rnames.keep, ]

# Process downstream ------------------------------------------------------

m.act.col <- ggplot(dat.louv.lst[[jmark1]], aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark1)
print(m.act.col)

m.act.numb <- ggplot(dat.louv.lst[[jmark1]], aes(x = umap1, y = umap2, label = cluster)) + geom_text() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark1)
print(m.act.numb)

m.repress.numb <- ggplot(dat.louv.lst[[jmark2]], aes(x = umap1, y = umap2, label = cluster)) + geom_text() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark2)
print(m.repress.numb)

m.repress.col <- ggplot(dat.louv.lst[[jmark2]], aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark2)
print(m.repress.col)

if (endsWith(inf.dbl.lda, ".rds")){
  out.lda <- readRDS(inf.dbl.lda)  # out.lda
} else {
  load(inf.dbl.lda, v=T)
}

tm.result <- posterior(out.lda)
# plot UMAP of double data
topics.mat <- tm.result$topics
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.dbl <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], mark = jmark.dbl, stringsAsFactors = FALSE)
dat.umap.dbl <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.dbl)

ggplot(dat.umap.dbl, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark.dbl.dash)


# Process fits ------------------------------------------------------------

# handle zeros optionally
w.lst <- sapply(act.repress.coord.lst, function(x){
  return(x$w)
})
plot(hist(w.lst))  # more active than repressive? Why is that? 

if (args$RemoveBadMixings){
  print("Number of cells in model before removing bad mixings:")
  print(length(act.repress.coord.lst))

  cells.remove.i <- which(w.lst <= 0.01 | w.lst >= 0.99)
  cells.remove <- names(act.repress.coord.lst)[cells.remove.i]
  act.repress.coord.removed.lst <- list()
  for (jcell in cells.remove){
    print(paste("Removing", jcell, "because w=", w.lst[[jcell]]))
    act.repress.coord.removed.lst[[jcell]] <- act.repress.coord.lst[[jcell]]
    act.repress.coord.lst[[jcell]] <- NULL
  }
  print("Number of cells in model after removing bad mixings:")
  print(length(act.repress.coord.lst))

  # remove cells from raw counts
  cells.keep <- names(act.repress.coord.lst)
  print("Dim of raw count mat before bad mixings:")
  print(dim(count.mat.dbl))
  count.mat.dbl <- count.mat.dbl[, cells.keep]
  print("Dim of raw count mat after bad mixings:")
  print(dim(count.mat.dbl))
}

# process fits
fits.out <- act.repress.coord.lst

out.dat <- data.frame(cell = names(w.lst), w = w.lst, stringsAsFactors = FALSE)
out.dat$experi <- sapply(as.character(out.dat$cell), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_"))

dat.umap.dbl.merge <- left_join(dat.umap.dbl, out.dat)

m.w <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = w)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

print(m.w)

# bimodal???
m.hist <- ggplot(out.dat, aes(x = w)) + geom_histogram() + 
  
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

print(m.hist)



# get active, repress index for each cell
# if clusters are now from clusters need eto rethink jcoord
cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  
  # rows are active, columns are repress I THINK?
  # TODO: assumes underscores be careful!
  jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
  jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
  
  if (grepl("_", jlouv.act)){
    jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
  } 
  if (grepl("_", jlouv.repress)){
    jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
  } 
  
  out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()


# set coords.dbl as factor: kind of old but keep it for historical reasons 
coords.dbl$louv.act.col <- sapply(as.factor(coords.dbl$louv.act), function(x) cbPalette[[x]])
coords.dbl$louv.repress.col <- sapply(as.factor(coords.dbl$louv.repress), function(x) cbPalette[[x]])

dat.umap.dbl <- left_join(dat.umap.dbl, coords.dbl)

dat.umap.dbl$louv.repress.char <- dat.umap.dbl$louv.repress
dat.umap.dbl$louv.active.char <- dat.umap.dbl$louv.act
dat.umap.dbl$pair <- paste(dat.umap.dbl$louv.act, dat.umap.dbl$louv.repress, sep = "_")

m.dbl.infer.colrepr.numbact <- ggplot(dat.umap.dbl, aes(x = umap1, y = umap2, color = louv.repress.col, label = louv.act)) + 
  geom_text(size = jtxtsize) + 
  scale_color_identity() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Colors are repressed clusters, numbers are active clusters")
# print(m.dbl.infer.colrepr.numbact)

m.dbl.infer.numbrepr.colact <- ggplot(dat.umap.dbl, aes(x = umap1, y = umap2, color = louv.act.col, label = louv.repress)) + 
  geom_text(size = jtxtsize) + 
  scale_color_identity() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("Colors are active clusters, numbers are repressed clusters") 
# print(m.dbl.infer.numbrepr.colact)

# fraction of cells in double marks

N <- nrow(dat.umap.dbl)
dat.umap.dbl.sum <- dat.umap.dbl %>%
  group_by(louv.act) %>%
  summarise(frac.cells = length(louv.act) / N) 

dat.umap.dbl.sum.repress <- dat.umap.dbl %>%
  group_by(louv.repress) %>%
  summarise(frac.cells = length(louv.repress) / N) %>%
  ungroup()


m.histo <- ggplot(coords.dbl, aes(x = w)) + geom_histogram(bins = 100) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(outpdf, useDingbats = FALSE)
ggplot(dat.umap.dbl.sum, aes(x = louv.act, y = frac.cells)) + geom_bar(stat = "identity") + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.umap.dbl.sum.repress, aes(x = louv.repress, y = frac.cells)) + geom_bar(stat = "identity") + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.histo)
lapply(list(m.act.col, m.repress.numb, m.dbl.infer.numbrepr.colact), print)
lapply(list(m.act.numb, m.repress.col, m.dbl.infer.colrepr.numbact), print)
multiplot(m.act.col + theme(legend.position = "none"), m.dbl.infer.numbrepr.colact, cols = 2)
multiplot(m.repress.col + theme(legend.position = "none"), m.dbl.infer.colrepr.numbact, cols = 2)
print(m.w)
print(m.hist)
print(m.act.col)
print(m.act.numb)
print(m.repress.col)
print(m.repress.numb)
print(m.dbl.infer.colrepr.numbact)
print(m.dbl.infer.numbrepr.colact)
dev.off()


# Get pairs ---------------------------------------------------------------

coords.dbl.stringent <- subset(coords.dbl, lnprob > logpcutoff)

# count all pairs
dbl.pairs <- paste(coords.dbl.stringent$louv.act, coords.dbl.stringent$louv.repress, sep = "_")
dbl.pairs.counts <- sort(table(dbl.pairs), decreasing = TRUE)
dbl.pairs.counts.filt <- dbl.pairs.counts[which(dbl.pairs.counts >= 10)]

plot(hist(dbl.pairs.counts, breaks = 100))



# Plot pairs --------------------------------------------------------------


pdf(outpdf.pairs, useDingbats = FALSE)

for (i in seq(length(dbl.pairs.counts))){
# for (i in seq(length(dbl.pairs.counts.filt))){
  
  jpair <- names(dbl.pairs.counts)[[i]]
  jcounts <- dbl.pairs.counts[[i]]
  print(paste(jpair, jcounts))
  # TODO: assumes underscores be careful!
  clst1 <- strsplit(jpair, split = "_")[[1]][[1]]
  clst2 <- strsplit(jpair, split = "_")[[1]][[2]]
  
  m.dbl <- ggplot(dat.umap.dbl %>% mutate(is.pair = pair == jpair), aes(x = umap1, y = umap2, color = is.pair)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark.dbl.dash, "double mark cluster pairs:", jpair, "N cells:", jcounts)) + 
    scale_color_manual(values = cbPalette)
  m.act.col.pair <- ggplot(dat.louv.lst[[jmark1]] %>% mutate(is.pair = cluster == clst1), 
                           aes(x = umap1, y = umap2, color = is.pair)) + geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jmark1, "Pair:", jpair))
  m.repress.col.pair <- ggplot(dat.louv.lst[[jmark2]] %>% mutate(is.pair = cluster == clst2), 
                               aes(x = umap1, y = umap2, color = is.pair)) + geom_point() + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    ggtitle(jmark2)
  # plot K36me3 and K27me3 side by side 
  print(m.dbl)
  multiplot(m.act.col.pair, m.repress.col.pair, cols = 2)
}
dev.off()



# Visualize matrix of specific pairs?  ------------------------------------

pdf(outpdf.mats, useDingbats = FALSE)
# for topic27, what are the repress topics?
jclsts <- unique(sapply(names(dbl.pairs.counts), function(x) strsplit(x, "_")[[1]][[1]]))
for (jclst in jclsts){
  print(jclst)
  cells.sub <- subset(coords.dbl, louv.act == jclst)$cell
  ll.mats <- lapply(cells.sub, function(jcell){
    jfit <- fits.out[[jcell]]$ll.mat
    return(jfit)
  })
  
  ll.mats.sum <- Reduce("+", ll.mats) / length(ll.mats)
  p.mats.sum <- SoftMax(ll.mats.sum)
  
  # use ggplot2
  ll.mats.sum.long <- tidyr::gather(data.frame(active = rownames(ll.mats.sum), ll.mats.sum), key = "repress", value = "minusLL", -active)
  p.mats.sum.long <- tidyr::gather(data.frame(active = rownames(ll.mats.sum), p.mats.sum), key = "repress", value = "logProb", -active)
  
  m.mats <- ggplot(ll.mats.sum.long, aes(x = active, y = repress, fill = minusLL)) + geom_tile() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_viridis_c(direction = 1) + ggtitle(paste(jmark1, jclst, collapse = " "))
  print(m.mats)
  m.pmats <- ggplot(p.mats.sum.long, aes(x = active, y = repress, fill = exp(logProb))) + geom_tile() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_viridis_c(direction = 1) + ggtitle(paste(jmark1, jclst, collapse = " "))
  print(m.pmats)
}
dev.off()



# Optinally add annotations -------------------------------------------------------------

# annotes
# can add additional annotations here

# dat.umap.dbl.annot <- left_join(dat.umap.dbl, annots.dat)
# dat.umap.dbl.annot <- dat.umap.dbl  # for unfixed no annots.dat

if (any(is.na(dat.umap.dbl$w))){
  print("There are double cells that are not assigned, likely because you removed NAs in the double data, removing these NAs for further analyses...")
  print(paste("Number of cells before filtering out NAs in double mat:", nrow(dat.umap.dbl)))
  print(head(dat.umap.dbl))
  dat.umap.dbl.annot <- subset(dat.umap.dbl, !is.na(w))  # for unfixed no annots.dat
  print(paste("Number of cells after filtering out NAs in double mat:", nrow(dat.umap.dbl.annot)))
  print(head(dat.umap.dbl.annot))
} else {
  print("No NAs detected in double marks, keeping all cells")
  dat.umap.dbl.annot <-dat.umap.dbl
  print(paste("Nubmer of cells in double mat:", nrow(dat.umap.dbl.annot)))
}

# factor as alphabetical?
# dat.umap.dbl.annot$celltype <- factor(dat.umap.dbl.annot$celltype, levels = sort(unique(dat.umap.dbl.annot$celltype)))

# dat.louv.lst.annot <- lapply(dat.louv.lst, function(x){
#   xdat <- left_join(x, annots.dat)
#   xdat$celltype <- factor(xdat$celltype, levels = sort(unique(xdat$celltype)))
#   return(xdat)
# })
dat.louv.lst.annot <- dat.louv.lst


# remove NAs
# dat.louv.lst.annot <- subset(dat.louv.lst.annot, !is.na(w))



# Unmix for all cells  ----------------------------------------------------


all.cells <- dat.umap.dbl.annot$cell
names(all.cells) <- all.cells
# lapply(seq_len(ncol(x)), function(i) x[,i])
col.i <- seq_len(ncol(count.mat.dbl))
names(col.i) <- colnames(count.mat.dbl)
all.x.raw <- lapply(col.i, function(i) count.mat.dbl[, i])  # https://stackoverflow.com/questions/6819804/how-to-convert-a-matrix-to-a-list-of-column-vectors-in-r/6823557

all.mixweights <- dat.umap.dbl.annot$w
names(all.mixweights) <- all.cells
all.louv.active <- dat.umap.dbl.annot$louv.active.char
all.louv.repress <- dat.umap.dbl.annot$louv.repress.char
names(all.louv.active) <- all.cells
names(all.louv.repress) <- all.cells
all.p.active <- lapply(all.louv.active, function(clstr.active) dat.impute1[clstr.active, ])
all.p.repress <- lapply(all.louv.repress, function(clstr.repress) dat.impute2.lst[[clstr.repress]])


jnames.all <- lapply(X = list(all.x.raw, all.mixweights, all.p.active, all.p.repress), FUN = names)
assertthat::assert_that(all(sapply(jnames.all, identical, jnames.all[[1]])))
# assertthat::assert_that(all(names(all.x.raw) == names(all.mixweights)))

system.time(
  x.raw.unmixed <- lapply(all.cells, function(jcell){
    # print(jcell)
    return(UnmixRawCounts(x.raw = all.x.raw[[jcell]], mixweight = all.mixweights[[jcell]], p.active = all.p.active[[jcell]], p.repress = all.p.repress[[jcell]], random.seed = 0))
  })
)
save(x.raw.unmixed, file = outf.unmixing)

hist(x.raw.unmixed[[1]]$p.cell.active.weights, col = 'red', main = paste0("Unmixing ", jmark.dbl, " to ", jmark1, "\n", all.cells[[1]]), xlab = paste0("Probability of Bin Assigned to ", jmark1))
hist(x.raw.unmixed[[2]]$p.cell.active.weights, col = 'red', main = paste0("Unmixing ", jmark.dbl, " to ", jmark1, "\n", all.cells[[2]]), xlab = paste0("Probability of Bin Assigned to ", jmark1))

# genomewide?
p.all <- unlist(lapply(x.raw.unmixed, function(sublst) sublst$p.cell.active.weights))
# qplot(x = p.all, geom = "density") + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle("Distribution of Probabilities Across All Cells")
plot(density(sample(x = p.all, size = 0.1 * length(p.all), replace = FALSE)), main = "Distribution of Probabilities Across All Cells", col = 'blue')
# hist(p.all, col = "red")


# create matrix of probabilities for MF
mat.probs <- as.data.frame(lapply(x.raw.unmixed, function(sublst) sublst$p.cell.active.weights))
# fix column names 
colnames(mat.probs) <- gsub("\\.", "-", colnames(mat.probs))
fwrite(mat.probs, file = outf.unmixing.probs.mat, sep = "\t", row.names = TRUE)

# jcmd <- paste("gzip", outf.unmixing.probs.mat)
# system(command = jcmd)

# View downstream ---------------------------------------------------------

pdf(outpdf.split, useDingbats = FALSE)

jcell <- names(all.x.raw)[[1]]
par(mfrow=c(3,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(all.x.raw[[1]][1:2000], type = "l")
plot(x.raw.unmixed[[1]]$x.raw.active[1:2000], type = "l")
plot(x.raw.unmixed[[1]]$x.raw.repress[1:2000], type = "l")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

x.sum <- all.x.raw[[jcell]]
x.act <- x.raw.unmixed[[jcell]]$x.raw.active
x.rep <- x.raw.unmixed[[jcell]]$x.raw.repress

plot(x.act + x.rep, x.sum)

# do umap to show this is legit
rnames <- rownames(count.mat.dbl)
x.act.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.active)), row.names = rnames)
x.repress.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.repress)), row.names = rnames)



count.mat.act.withunmixed <- cbind(count.mat.act, as.matrix(x.act.mat))
# check colnames
unique(sapply(colnames(count.mat.act.withunmixed), function(x) ClipLast(x, jsep = "_")))
# fix cnames, dots to dashes
colnames(count.mat.act.withunmixed) <- gsub("\\.", "-", colnames(count.mat.act.withunmixed))

count.mat.repress.withunmixed <- cbind(count.mat.repress, as.matrix(x.repress.mat))
# check
unique(sapply(colnames(count.mat.repress.withunmixed), function(x) ClipLast(x, jsep = "_")))  # no need to fix 
colnames(count.mat.repress.withunmixed) <- gsub("\\.", "-", colnames(count.mat.repress.withunmixed))

# get clean mat for output
x.active.mat.clean <- as.matrix(x.act.mat)
colnames(x.active.mat.clean) <- gsub("\\.", "-", colnames(x.active.mat.clean))
x.repress.mat.clean <- as.matrix(x.repress.mat)
colnames(x.repress.mat.clean) <- gsub("\\.", "-", colnames(x.repress.mat.clean))



lsi.act.out <- RunLSI(as.matrix(count.mat.act.withunmixed))
# lsi.act.out <- RunLSI(as.matrix(count.mat.act))

umap.out.lsi.act <- umap(lsi.act.out$u, config = jsettings)
dat.umap.long.lsi.act <- data.frame(cell = rownames(umap.out.lsi.act[["layout"]]), umap1 = umap.out.lsi.act[["layout"]][, 1], umap2 = umap.out.lsi.act[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long.lsi.act <- DoLouvain(lsi.act.out$u, custom.settings.louv = jsettings, dat.umap.long.lsi.act)

dat.umap.long.lsi.act.merge <- dat.umap.long.lsi.act %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))
# left_join(., annots.dat)

ggplot(dat.umap.long.lsi.act.merge, aes(x = umap1, y = umap2, color = louvain, shape = plate)) + geom_point(size = 2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 

ggplot(dat.umap.long.lsi.act.merge, aes(x = umap1, y = umap2, color = louvain, shape = plate)) + geom_point(size = 2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)  + facet_wrap(~plate)


# do other mark
lsi.repress.out <- RunLSI(as.matrix(count.mat.repress.withunmixed))

umap.out.lsi.repress <- umap(lsi.repress.out$u, config = jsettings)
dat.umap.long.lsi.repress <- data.frame(cell = rownames(umap.out.lsi.repress[["layout"]]), umap1 = umap.out.lsi.repress[["layout"]][, 1], umap2 = umap.out.lsi.repress[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long.lsi.repress <- DoLouvain(lsi.repress.out$u, custom.settings.louv = jsettings, dat.umap.long.lsi.repress)

dat.umap.long.lsi.repress.merge <- dat.umap.long.lsi.repress %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))
# left_join(., annots.dat)

ggplot(dat.umap.long.lsi.repress.merge, aes(x = umap1, y = umap2, color = louvain, shape = plate)) + geom_point(size = 2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) 

ggplot(dat.umap.long.lsi.repress.merge, aes(x = umap1, y = umap2, color = louvain, shape = plate)) + geom_point(size = 2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)  + facet_wrap(~plate)


cells.random <- sample(all.cells, size = 10)
for (cell in cells.random){
  hist(x.raw.unmixed[[cell]]$p.cell.active.weights, col = 'red', main = paste0("Unmixing ", jmark.dbl, " to ", jmark1, "\n", cells.random[[cell]]), xlab = paste0("Probability of Bin Assigned to ", jmark1))
}

dev.off()

# Write matrix to output for LDA  -----------------------------------------

# not sure if works for LDA if rows all empty?  let LDA handle that
saveRDS(count.mat.act.withunmixed, file = outf.mergedmat.act)
saveRDS(count.mat.repress.withunmixed, file = outf.mergedmat.repress)

# write unmixed only
# 2020-10-22 need to subset to be same rownames as respectivve LDAs
print("Writing unmixed only, subsetting to same rownames as respective LDAs for projecting")
print(head(rnames.act))
print(head(rownames(x.active.mat.clean)))

print(all(rownames(x.active.mat.clean) %in% rnames.act))
print(all(rnames.act %in% rownames(x.active.mat.clean)))

print(all(rownames(x.repress.mat.clean) %in% rnames.repress))
print(all(rnames.repress %in% rownames(x.repress.mat.clean)))

x.active.mat.clean.ForProj <- x.active.mat.clean[rnames.act, ]
x.repress.mat.clean.ForProj <- x.repress.mat.clean[rnames.repress, ]

# check dims before and after
print("Dim of mark1 before:")
print(dim(x.active.mat.clean))
print("Dim of mark1 after:")
print(dim(x.active.mat.clean.ForProj))

print("Dim of mark2 before:")
print(dim(x.repress.mat.clean))
print("Dim of mark2 after:")
print(dim(x.repress.mat.clean.ForProj))

saveRDS(x.active.mat.clean.ForProj, file = outf.unmixedmat.act)
saveRDS(x.repress.mat.clean.ForProj, file = outf.unmixedmat.repress)

print(x.active.mat.clean[1:5, 1:5])
print(x.repress.mat.clean[1:5, 1:5])

# Write bad cells to file
if (args$RemoveBadMixings){
  fwrite(data.frame(cell = names(w.lst)[cells.remove.i], w = w.lst[cells.remove.i]), file = outtxt.badcells)
  save(cells.remove, act.repress.coord.removed.lst, file = outrdata.badcells)
}


