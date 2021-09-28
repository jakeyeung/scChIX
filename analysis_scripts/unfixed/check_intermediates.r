# Jake Yeung
# Date of Creation: 2021-09-13
# File: ~/projects/scChIX/analysis_scripts/unfixed/check_intermediates.r
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(irlba)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Set variables from args -------------------------------------------------

args <- list()

jmark1 <- "K4m1"
jmark2 <- "K27m3"
jmarks.pair <- c(jmark1, jmark2)
jmark.dbl.dash <- paste(jmarks.pair, collapse = "-")
jmark.dbl <- paste(c(jmark1, jmark2), collapse = "-")
names(jmarks.pair) <- jmarks.pair

inf.in <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/scchix_inputs_objs/scchix_inputs_clstr_by_celltype_K4m1-K27m3.removeNA_FALSE.RData")
inf.out <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_K4m1-K27m3.RData")
inf.act <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/countmat_output_filt.K4m1.rds")
inf.rep <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/countmat_output_filt.K27m3.rds")
inf.dbl.lda <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/LDA_outputs_init/ldaOut.countmat_var_filt.K4m1-K27m3.Robj")





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




load(inf.dbl.lda, v=T)
tm.result <- posterior(out.lda)
# plot UMAP of double data
topics.mat <- tm.result$topics
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.dbl <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], mark = jmark.dbl, stringsAsFactors = FALSE)
dat.umap.dbl <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.dbl)








# Process fits ------------------------------------------------------------

# handle zeros optionally
w.lst <- sapply(act.repress.coord.lst, function(x){
  return(x$w)
})
plot(hist(w.lst))  # more active than repressive? Why is that?



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







jtxtsize <- 2
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


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
print(m.dbl.infer.colrepr.numbact)

m.dbl.infer.numbrepr.colact <- ggplot(dat.umap.dbl, aes(x = umap1, y = umap2, color = louv.act.col, label = louv.repress)) +
  geom_text(size = jtxtsize) +
  scale_color_identity() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  ggtitle("Colors are active clusters, numbers are repressed clusters")
print(m.dbl.infer.numbrepr.colact)

# fraction of cells in double marks

N <- nrow(dat.umap.dbl)
dat.umap.dbl.sum <- dat.umap.dbl %>%
  group_by(louv.act) %>%
  summarise(frac.cells = length(louv.act) / N)

dat.umap.dbl.sum.repress <- dat.umap.dbl %>%
  group_by(louv.repress) %>%
  summarise(frac.cells = length(louv.repress) / N) %>%
  ungroup()




# Get pairs ---------------------------------------------------------------

pcutoff <- 0.95
logpcutoff <- log(pcutoff)
coords.dbl.stringent <- subset(coords.dbl, lnprob > logpcutoff)

# count all pairs
dbl.pairs <- paste(coords.dbl.stringent$louv.act, coords.dbl.stringent$louv.repress, sep = "_")
dbl.pairs.counts <- sort(table(dbl.pairs), decreasing = TRUE)
dbl.pairs.counts.filt <- dbl.pairs.counts[which(dbl.pairs.counts >= 10)]

plot(hist(dbl.pairs.counts, breaks = 100))





# Unmix for all cells  ----------------------------------------------------

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



# create matrix of probabilities for MF
mat.probs <- as.data.frame(lapply(x.raw.unmixed, function(sublst) sublst$p.cell.active.weights))
# fix column names
colnames(mat.probs) <- gsub("\\.", "-", colnames(mat.probs))





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


inf.check <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_inputs/countmats/countmat_var_filt.K4m1-K27m3.rds")
mat.check <- readRDS(inf.check)

inf.ldacheck <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/LDA_outputs_init/ldaOut.countmat_var_filt.K4m1-K27m3.Robj")
load(inf.ldacheck, v=T)

inf.ldacheck2 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/lda_output_filt.K4m1-K27m3.rds")
ldacheck2 <- readRDS(inf.ldacheck2)

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


# Checkmore ---------------------------------------------------------------

inf.check1 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/lda_output_filt.K4m1.rds")
inf.check2 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/lda_output_filt.K27m3.rds")
inf.check3 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/objs_from_LDA/countmat_output_filt.K4m1-K27m3.rds")

assertthat::assert_that(file.exists(inf.check1))
assertthat::assert_that(file.exists(inf.check2))
assertthat::assert_that(file.exists(inf.check3))


check1 <- readRDS(inf.check1)
check2 <- readRDS(inf.check2)
check3 <- readRDS(inf.check3)

tm.result.lst <- list()
tm.result.lst[[jmark1]] <- posterior(check1)
tm.result.lst[[jmark2]] <- posterior(check2)


rnames.lst <- list(colnames(tm.result.lst[[jmark1]]$terms), colnames(tm.result.lst[[jmark2]]$terms), rownames(check3))
rnames.common <- Reduce(f = intersect, x = rnames.lst)

jcheck <- rnames.lst[[1]][!rnames.lst[[1]] %in% rnames.common]
