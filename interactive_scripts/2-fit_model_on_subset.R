# Jake Yeung
# Date of Creation: 2021-04-15
# File: ~/projects/scChIX/interactive_scripts/2-fit_-model_on_subset.R
#

rm(list=ls())

library(scChIX)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(igraph)
library(umap)


timestart <- Sys.time()

ncores <- 1

# bounds for fitting mixing fraction
wlower <- 0
wupper <- 1

# method for fitting mixing fraction
jmethod <- "Brent"

# load objs

data(ObjsForFittingSubset)

count.mat.dbl <- count.mat.dbl.subset

cell.count.raw.merged.lst <- as.list(as.data.frame(as.matrix(count.mat.dbl)))

system.time(
  # act.repress.coord.lst <- parallel::mclapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
  act.repress.coord.lst <- lapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
    optim.out <- FitMixingWeight(cell.count.raw.merged = cell.count.raw.merged,
                                 dat.impute.repress.lst = dat.impute.repress.lst,
                                 dat.impute.active = dat.impute.active, w.init = 0.5, w.lower = wlower, w.upper = wupper, jmethod = "Brent")
    ll.mat <- GetLLMerged(optim.out$par, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = TRUE)
    return(list(ll.mat = ll.mat, w = optim.out$par))
  # }, mc.cores = ncores)
  })
)



# Visualize a cell  -------------------------------------------------------

# data(GroundTruthLabels)

annots.dat <- LoadCellAnnotsEtOH.MorePlates() %>%
  dplyr::rename(cell.makenames = cell,
                cell = cell.orig)

(jcell.tmp <- names(act.repress.coord.lst)[[1]])
(jclst <- subset(annots.dat, cell == jcell.tmp)$celltype)

jcell.tmp <- "BM-EtOH-MNctrl-K27m3-K9m3-p1_220"
jclst <- "Granulocyte"  # ground truth

LL.mat.tmp <- act.repress.coord.lst[[jcell.tmp]]$ll.mat
LL.mat.tmp <- data.frame(H3K27me3 = rownames(LL.mat.tmp), LL.mat.tmp, stringsAsFactors = FALSE)
LL.long.tmp <- data.table::melt(LL.mat.tmp, id.vars = "H3K27me3", value.name = "LL", variable.name = "H3K9me3")

k27me3.ordering.all <- c("topic16", "topic3", "topic9", "topic25", "topic11")
k9me3.ordering.all <- c("topic12", "topic28", "topic1", "topic20")
ctype.ordering.all <- c("Granulocytes", "Bcells", "NKcells")

k27me3.clstr.labs <- list("topic16" = "Granulocytes",
                          "topic3" = "Bcells",
                          "topic9"  = "Bcells",
                          "topic25" = "NKcells",
                          "topic11" = "NKcells")
k9me3.clstr.labs <- list("topic12" = "Granulocytes",
                          "topic28" = "Bcells",
                          "topic1"  = "NKcells",
                          "topic20" = "NKcells")



LL.long.tmp$H3K27me3 <- factor(LL.long.tmp$H3K27me, levels = k27me3.ordering.all)
LL.long.tmp$H3K9me3 <- factor(LL.long.tmp$H3K9me3, levels = k9me3.ordering.all)

m <- ggplot(LL.long.tmp, aes(x = H3K27me3, y = H3K9me3, fill = LL)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c() +
  ggtitle(paste("Ground truth:", jclst, "\nCell name: ", jcell.tmp)) +
  theme(aspect.ratio= 4/8, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom")
print(m)



# Summarize all cells  ----------------------------------------------------


# if louvains are now from clusters need eto rethink jcoord
cell.vec <- names(act.repress.coord.lst)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- act.repress.coord.lst[[jcell]]
  jweight <- act.repress.coord.lst[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)

  jclstr.k27me3 <- rownames(p.mat)[[jcoord[[1]]]]
  jclstr.k9me3 <- colnames(p.mat)[[jcoord[[2]]]]

  if (grepl("_", jclstr.k27me3)){
    jclstr.k27me3 <- strsplit(jclstr.k27me3, split = "_")[[1]][[2]]
  }
  if (grepl("_", jclstr.k9me3)){
    jclstr.k9me3 <- strsplit(jclstr.k9me3, split = "_")[[1]][[2]]
  }

  out.dat <- data.frame(cell = jcell, clstr.k27me3 = jclstr.k27me3, clstr.k9me3 = jclstr.k9me3, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()


coords.dbl.annots <- left_join(coords.dbl, annots.dat)
coords.dbl.annots$clstr.k27me3 <- factor(coords.dbl.annots$clstr.k27me3, levels = k27me3.ordering.all)
coords.dbl.annots$clstr.k9me3 <- factor(coords.dbl.annots$clstr.k9me3, levels = k9me3.ordering.all)

# annotate clusters
coords.dbl.annots$clstr.k27me3.annot <- sapply(as.character(coords.dbl.annots$clstr.k27me3), function(clst){
  k27me3.clstr.labs[[clst]]
})
coords.dbl.annots$clstr.k9me3.annot <- sapply(as.character(coords.dbl.annots$clstr.k9me3), function(clst){
  k9me3.clstr.labs[[clst]]
})

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.grid <- ggplot(coords.dbl.annots, aes(x = clstr.k27me3.annot, y = clstr.k9me3.annot, color = celltype)) +
  geom_point(position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab("H3K27me3 Cell Type Clusters") + ylab("H3K9me3 Cell Type Clusters") +
  ggtitle("Each dot is a double stained cell, colored by ground truth label.\nX-Y shows the cluster pair it is assigned")
print(m.grid)


# Deconvolve double-incubated count matrix ------------------------------------------------------


all.cells <- coords.dbl$cell
names(all.cells) <- all.cells
col.i <- seq_len(ncol(count.mat.dbl))
names(col.i) <- colnames(count.mat.dbl)
all.x.raw <- lapply(col.i, function(i) count.mat.dbl[, i])  # https://stackoverflow.com/questions/6819804/how-to-convert-a-matrix-to-a-list-of-column-vectors-in-r/6823557

all.mixweights <- coords.dbl$w
names(all.mixweights) <- all.cells
all.clstr.k27me3 <- as.character(coords.dbl$clstr.k27me3)
all.clstr.k9me3 <- as.character(coords.dbl$clstr.k9me3)
names(all.clstr.k27me3) <- all.cells
names(all.clstr.k9me3) <- all.cells
all.p.active <- lapply(all.clstr.k27me3, function(clstr.k27me3) dat.impute.active[clstr.k27me3, ])
all.p.repress <- lapply(all.clstr.k9me3, function(clstr.k9me3) dat.impute.repress.lst[[clstr.k9me3]])


jnames.all <- lapply(X = list(all.x.raw, all.mixweights, all.p.active, all.p.repress), FUN = names)
assertthat::assert_that(all(sapply(jnames.all, identical, jnames.all[[1]])))
# assertthat::assert_that(all(names(all.x.raw) == names(all.mixweights)))

system.time(
  x.raw.unmixed <- lapply(all.cells, function(jcell){
    # print(jcell)
    return(UnmixRawCounts(x.raw = all.x.raw[[jcell]], mixweight = all.mixweights[[jcell]], p.active = all.p.active[[jcell]], p.repress = all.p.repress[[jcell]], random.seed = 0))
  })
)


# Project onto UMAP  ------------------------------------------------------


# load LDAs
# load("/home/jyeung/projects/scChIX/data/H3K27me3_LdaOutputs.RData", v=T)
# load("/home/jyeung/projects/scChIX/data/H3K9me3_LdaOutputs.RData", v=T)

data(H3K27me3_LdaOutputs)  # out.lda.k27me3
data(H3K9me3_LdaOutputs)  # out.lda.k9me3

rnames <- rownames(count.mat.dbl)
x.k27me3.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.active)), row.names = rnames)
x.k9me3.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.repress)), row.names = rnames)

set.seed(0)
nsubset <- 10  # faster
rand.i1 <- sample(seq(ncol(x.k27me3.mat)), size = nsubset)
rand.i2 <- sample(seq(ncol(x.k9me3.mat)), size = nsubset)

print("Projecting to a few unmixed cells to K27me3 manifold")
system.time(
  out.lda.predict.k27me3 <- topicmodels::posterior(out.lda.k27me3, t(as.matrix(x.k27me3.mat[, rand.i1])))
)

print("Projecting to a few unmixed cells to K9me3 manifold")
system.time(
  out.lda.predict.k9me3 <- topicmodels::posterior(out.lda.k27me3, t(as.matrix(x.k27me3.mat[, rand.i2])))
)


# Plot outputs -----------------------------------------------------------

# jsuffix <- ".FewerTopics2"
# jprefix.base <- "TopBins_autosomesOnly"
# jprefix <- paste0(jprefix.base, jsuffix, "_KeepAllPlates_clstrbytopics")
# jprefix.input <- jprefix
# jprefix.output <- paste0(jprefix, "_0_1")
# jprefix.ldaoutput <- paste0(jprefix, "_0_1_RemoveBadMixings")
#
# jmarks <- c("K27m3", "K9m3")
# names(jmarks) <- jmarks
# indir.lda <-  paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/projections_LDA_outputs/afterUnmixing_", jprefix.ldaoutput, "_K27m-K9m3")
#

jsettings <- umap::umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("K27m3", "K9m3")
names(jmarks) <- jmarks

# load("/home/jyeung/projects/scChIX/data/H3K27me3_ProjectionOutput.RData", v=T)
# load("/home/jyeung/projects/scChIX/data/H3K9me3_ProjectionOutput.RData", v=T)

data(H3K27me3_ProjectionOutput.RData)  # out.proj.K27m3
data(H3K9me3_ProjectionOutput.RData)  # out.proj.K9m3

out.proj.lst <- list("K27m3" = out.proj.K27m3,
                     "K9m3" = out.proj.K9m3)

dat.umap.merge <- lapply(jmarks, function(jmark){
  out.proj.mark <- out.proj.lst[[jmark]]
  # inf.lda <- file.path(indir.lda, paste0("project_unmixed_", jmark, ".RData"))
  # load(inf.lda)
  #  tm.result <- topicmodels::posterior(out.objs$out.lda)
  tm.result <- topicmodels::posterior(out.proj.mark$out.lda)

  topics.mat.orig <- tm.result$topics
  topics.mat.proj <- out.proj.mark$out.lda.predict$topics

  umap.out <- umap::umap(topics.mat.orig, config = jsettings)
  rownames(umap.out$layout) <- rownames(topics.mat.orig)
  dat.umap <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  umap.proj <- predict(umap.out, topics.mat.proj)
  dat.umap.proj <- data.frame(cell = rownames(umap.proj), umap1 = umap.proj[, 1], umap2 = umap.proj[, 2], stringsAsFactors = FALSE)

  dat.umap.merge <- rbind(dat.umap %>% mutate(cond = "single"), dat.umap.proj %>% mutate(cond = "double")) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_", jsep.out = "_")) %>%
    left_join(., annots.dat)
  dat.umap.merge$mark <- jmark
  return(dat.umap.merge)
})


# merge the single umaps into one, seaprate by umap2
dat.merge.singles <- dat.umap.merge %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(umap1 = scale(umap1, center = TRUE, scale = TRUE),
         umap2 = scale(umap2, center = TRUE, scale = TRUE))

yrange <- diff(range(dat.merge.singles$umap2)) / 1.5
xrange <- diff(range(dat.merge.singles$umap1)) / 1.5

dat.merge.singles <- dat.merge.singles %>%
  rowwise() %>%
  mutate(umap2.shift = umap2,
         umap1.shift = ifelse(mark == jmarks[[1]], umap1 + xrange, umap1 - xrange))

m.umaps.dbls.ctype <- lapply(jmarks, function(jmark){
  jsub <- subset(dat.merge.singles, mark == jmark)
  m <- ggplot(jsub, aes(x = umap1, y = umap2, color = celltype, shape = cond)) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, "Single + Dbl Split")) +
    facet_wrap(~cond)
  return(m)
})
print(m.umaps.dbls.ctype)

m <- ggplot(dat.merge.singles, aes(x = -1 * umap1.shift, y = umap2.shift, color = celltype, shape = mark, group = cell)) +
  geom_point() +
  geom_path(size = 0.1, alpha = 0.1) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggtitle(paste( "Linked map of H3K27me3 (left) and H3K9me3 (right)")) +
  xlab("UMAP1 (shifted)") + ylab("UMAP2")
print(m)


print(Sys.time() - timestart)

