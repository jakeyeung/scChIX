# Jake Yeung
# Date of Creation: 2021-08-19
# File: ~/projects/scChIX/analysis_scripts/5-transfer_labels_from_scChIX.K9me3.check_K36.refilt.transfer_twice.again.R
#




rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(ggforce)

library(topicmodels)

library(scchicFuncs)

library(JFuncs)

jsize <- 1
jspread <- 7

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- jspread


hubprefix <- "/home/jyeung/hub_oudenaarden"


# Functions ---------------------------------------------------------------


CellToNearestNeighbors <- function(xvec, xmat, jmeth = "euclidean"){
  assertthat::assert_that(jmeth == "euclidean")
  xdist <- apply(xmat, 1, function(jrow) sum((xvec - jrow) ^ 2))
  xdist.sort <- sort(xdist, decreasing=FALSE)
  return(xdist.sort)
}

nnearest <- 30
nnearest2 <- 5

# Load matrix  ------------------------------------------------------------

jdate <- "2021-08-19"
jquant <- "manual2nocenterfilt2"


jmark1 <- "K36"; jmark2 <- "K9m3"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")

names(jmarks) <- jmarks
jmarkdbl <- paste(c(jmark1, jmark2), collapse = "-")

jstr <- paste(c(jmarks, jmarkdbl), collapse = "_")

jprefix <- "var_filtered"
jname <- paste(jprefix, jquant, jstr, sep = "_")

infrdata <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_outputs_objs/", jname, "/unmix_scchix_inputs_clstr_by_celltype_", jmarkdbl, ".removeNA_FALSE.RData"))
assertthat::assert_that(file.exists(infrdata))



# jsuffix <- "dbl_k36_k9m3_cleaned"
infs <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline/", jname, "/", jstr, "/Gastru_Unmixed_DblMark.", jname, ".", jmark, ".RData")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


# load mats
# metamain <- file.path(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_", jprefix))
metamain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA/", jname))
assertthat::assert_that(dir.exists(metamain))

dat.meta.lst <- lapply(jmarks, function(jmark){
  infmeta <- file.path(metamain, paste0("celltyping_output_filt.", jmark, ".", jdate, ".rds"))
  print(infmeta)
  dat.meta <- readRDS(infmeta)
  return(dat.meta)
})




# Load data  --------------------------------------------------------------


out.lst.lst <- lapply(jmarks, function(jmarktmp){
  load(infs[[jmarktmp]], v=T)  # out.objs, out.lda.predict
  out.lst <- list(out.objs = out.objs, out.lda.predict = out.lda.predict)
})

tm.lst <- lapply(out.lst.lst, function(x) posterior(x$out.objs$out.lda))

umap.objs.lst <- lapply(jmarks, function(jmarktmp){
  tm.orig <- tm.lst[[jmarktmp]]
  umap.out <- umap(tm.orig$topics, config = jsettings)
  return(umap.out)
})

dat.umap.merge.lst <- lapply(jmarks, function(jmarktmp){
  umap.out <- umap.objs.lst[[jmarktmp]]
  tm.orig <- tm.lst[[jmarktmp]]
  out.lda.predict <- out.lst.lst[[jmarktmp]]$out.lda.predict
  dat.umap.orig <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  dat.umap.orig <- DoLouvain(topics.mat = tm.orig$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.orig)
  dat.umap.orig.annot <- left_join(dat.umap.orig %>% mutate(type = "single") %>% dplyr::select(-louvain),
                                   dat.meta.lst[[jmarktmp]] %>% dplyr::select(c(cell, cluster)), by = "cell")

  # add projections
  umap.out.pred.layout <- predict(umap.out, data = out.lda.predict$topics)
  dat.umap.pred.annot <- data.frame(cell = rownames(umap.out.pred.layout), umap1 = umap.out.pred.layout[, 1], umap2 = umap.out.pred.layout[, 2], stringsAsFactors = FALSE) %>%
    mutate(type = "dbl",
           cluster = "na")

  dat.umap.merge <- rbind(dat.umap.orig.annot, dat.umap.pred.annot)
  dat.umap.merge$mark <- jmarktmp
  dat.umap.merge <- dat.umap.merge %>%
    rowwise() %>%
    mutate(stage = as.character(strsplit(cell, split = "-")[[1]][[1]]))
  return(dat.umap.merge)
})


# Plot linked UMAP  -------------------------------------------------------



dat.merge.rbind <- bind_rows(dat.umap.merge.lst) %>%
  group_by(mark) %>%
  mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
         umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
         umap1.shift = ifelse(mark == "K36", umap1.scale - 5, umap1.scale + 5)) %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]])



ggplot(dat.merge.rbind, aes(x = umap1.shift, y = 1 * umap2.scale, group = cell, color = stage)) +
  geom_point() +
  ggtitle("Double + single cells") +
  geom_path(alpha = 0.05) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# Label by UMAP location ------------------------------------------------

jsub <- subset(dat.merge.rbind, type == "dbl" & mark == jmarks[[2]] & umap2.scale > -1)
cell2lab <- hash::hash(jsub$cell, jsub$umap1.shift)

dat.merge.rbind$label.cont <- sapply(dat.merge.rbind$cell, function(x) AssignHash(x = x, jhash = cell2lab, null.fill = NA))

ggplot(dat.merge.rbind %>% filter(type == "dbl") %>% arrange(label.cont), aes(x = umap1.shift, y = umap2.scale, color = label.cont, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.05) +
  theme_bw() +
  ggtitle(paste(jmarks, collapse = "_")) +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge.rbind %>% filter(type == "dbl" & umap1.shift > -2.5) %>% arrange(label.cont), aes(x = umap1.shift, y = umap2.scale, color = label.cont, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.05) +
  theme_bw() +
  ggtitle(paste(jmarks, collapse = "_")) +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Do transfer but one cell at a time  -------------------------------------


jmark <- jmarks[[1]]

topics.mat1 <- tm.lst[[jmark]]$topics
topics.mat2 <- out.lst.lst[[jmark]]$out.lda.predict$topics

xvec.lst <- lapply(seq_len(nrow(topics.mat2)), function(i) topics.mat2[i, ])
names(xvec.lst) <- rownames(topics.mat2)

system.time(
  xdist.sort.lst <- lapply(xvec.lst, function(jrow) CellToNearestNeighbors(jrow, topics.mat1))
)


# get labels
cell2cluster <- hash::hash(dat.merge.rbind$cell, dat.merge.rbind$cluster)
xlab.sort.lst <- lapply(xdist.sort.lst, function(xvec) sapply(names(xvec), function(x) AssignHash(x, cell2cluster, null.fill = NA)))

# vote
jnames <- names(xlab.sort.lst)
names(jnames) <- jnames
xlab.sum.lst <- lapply(jnames, function(jname){
  xvec <- xlab.sort.lst[[jname]]
  data.frame(table(xvec[1:nnearest])) %>%
    mutate(rnk = rank(-Freq, ties.method = "random"),
           Frac = Freq / sum(Freq)) %>%
    arrange(rnk) %>%
    mutate(cell = jname) %>%
    dplyr::rename(cluster = Var1)
})

# check for ties and break them randomly
xlab.top <- lapply(xlab.sum.lst, function(x) x[1, ]) %>%
  bind_rows()
xlab.top$cell <- names(xlab.sum.lst)
table(xlab.top$rnk)

# Show K36 single and double  ---------------------------------------------

cell2cluster.sum <- hash::hash(xlab.top$cell, xlab.top$cluster)

dat.merge.rbind.dbl <- subset(dat.merge.rbind, type == "dbl") %>%
  rowwise() %>%
  mutate(cluster = AssignHash(cell, cell2cluster.sum, cluster))

dat.merge.rbind.single <- subset(dat.merge.rbind, type == "single")

dat.merge.rbind2 <- rbind(dat.merge.rbind.single,dat.merge.rbind.dbl)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
assertthat::assert_that(length(cbPalette) == length(unique(cbPalette)))

m.clsts <- ggplot(dat.merge.rbind2 %>% filter(mark == jmarks[[1]]), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  # facet_wrap(~type) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Annotate clusters  ------------------------------------------------------



print(m.clsts)



# Check clsts  ------------------------------------------------------------




# Annotate from before  ---------------------------------------------------

inf.meta.proj <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/celltyping_K36_first_try.2021-08-02.txt"
dat.meta.proj <- fread(inf.meta.proj)

cell2celltype.proj <- hash::hash(dat.meta.proj$cell, dat.meta.proj$celltype)

jcheck <- subset(dat.merge.rbind2, mark == "K36" & type == "single")
jcheck$celltype.check <- sapply(jcheck$cell, function(x) AssignHash(x, cell2celltype.proj, null.fill = x))

m.check <- ggplot(jcheck, aes(x = umap1, y = umap2, color = celltype.check)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette)

multiplot(m.clsts, m.check, cols = 2)



# Label cell types  -------------------------------------------------------



# topic10 = cluster6 = "erythroid"
clst.annots <- list()
clst.annots["cluster9"] <- "Erythroid"
clst.annots["cluster2"] <- "Endothelial"
clst.annots["cluster8"] <- "Epithelial"
clst.annots["cluster4"] <- "Neurons"
clst.annots["cluster1"] <- "ConnectiveTissueProg"
clst.annots["cluster5"] <- "WhiteBloodCells"
clst.annots["cluster10"] <- "NeuralTubeNeuralProgs"
clst.annots["cluster7"] <- "Stromal"
clst.annots["cluster3"] <- "SchwannCellPrecusor"

clst.annots["cluster6"] <- "Stromal"  # two clusters of stromal

clst.annots.hash <- hash(clst.annots)

jclst <- "cluster10"
ggplot(dat.merge.rbind2 %>% filter(mark == jmarks[[1]]) %>% mutate(cluster = ifelse(cluster == jclst, clst.annots[[jclst]], "NotInCluster")),
       aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  ggtitle(jclst, clst.annots[[jclst]]) +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  # facet_wrap(~type) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# reannotate
dat.merge.k36 <- subset(dat.merge.rbind2, mark == "K36") %>%
  rowwise() %>%
  mutate(celltype = AssignHash(cluster, clst.annots, null.fill = cluster))

m.annot <- ggplot(dat.merge.k36, aes(x = umap1, y = umap2, color = celltype),
       aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  # facet_wrap(~type) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.annot)

dat.merge.k36.pretty <- dat.merge.k36 %>%
  dplyr::rename(cluster.louvain = cluster,
                cluster = celltype)

# Write meta data ---------------------------------------------------------

outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix"
outdir <- file.path(outmain, jname)
dir.create(outdir)

outf <- file.path(outdir, paste0("celltyping_", jmark, ".", Sys.Date(), ".filt2.spread_", jspread, ".txt"))
outpdf <- file.path(outdir, paste0("celltyping_", jmark, ".", Sys.Date(), ".filt2.spread_", jspread, ".pdf"))
fwrite(x = dat.merge.k36, file = outf, sep = "\t")

pdf(outpdf, useDingbats = FALSE)
  print(m.clsts)
  print(m.annot)
dev.off()



# Transfer k9m3 dbl to k9m3  --------------------------------------------

# measure each k9m3 cell based on the distance from the K9 part of split cells
# use distance in latent space

jmark2 <- jmarks[[2]]
topics.mat.ref <- out.lst.lst[[jmark2]]$out.lda.predict$topics
topics.mat.query <- tm.lst[[jmark2]]$topics

xvec.query.lst <- lapply(seq_len(nrow(topics.mat.query)), function(i) topics.mat.query[i, ])
names(xvec.query.lst) <- rownames(topics.mat.query)

system.time(
  xdist.sort.query.lst <- lapply(xvec.query.lst, function(jrow) CellToNearestNeighbors(jrow, topics.mat.ref))
)


# get labels
# cell2cluster.after <- hash::hash(dat.merge.rbind$cell, dat.merge.rbind$cluster)
cell2cluster.query <- hash::hash(dat.merge.k36.pretty$cell, dat.merge.k36.pretty$cluster)
xlab.sort.query.lst <- lapply(xdist.sort.query.lst, function(xvec) sapply(names(xvec), function(x) AssignHash(x, cell2cluster.query, null.fill = NA)))


# vote
jnames <- names(xlab.sort.query.lst)
names(jnames) <- jnames
xlab.sum.query.lst <- lapply(jnames, function(jname){
  xvec <- xlab.sort.query.lst[[jname]]
  data.frame(table(xvec[1:nnearest2])) %>%
    mutate(rnk = rank(-Freq, ties.method = "random"),
           Frac = Freq / sum(Freq)) %>%
    arrange(rnk) %>%
    mutate(cell = jname) %>%
    dplyr::rename(cluster = Var1)
})

# check for ties and break them randomly
xlab.query.top <- lapply(xlab.sum.query.lst, function(x) x[1, ]) %>%
  bind_rows()
xlab.query.top$cell <- names(xlab.sum.query.lst)
table(xlab.query.top$rnk)


# Check output ------------------------------------------------------------

dat.merge.k9m3.single <- subset(dat.merge.rbind2, mark == jmark2 & type == "single", select = c(-cluster)) %>%
  left_join(., subset(xlab.query.top, select = c(cell, cluster)))

dat.merge.k9m3.dbl <- subset(dat.merge.rbind2, mark == jmark2 & type == "dbl", select = c(-cluster)) %>%
  left_join(., subset(dat.merge.k36.pretty, select = c(cell, cluster)))

# rename cluster
dat.merge.k9m3.merged <- rbind(dat.merge.k9m3.single, dat.merge.k9m3.dbl)

m.clsts2 <- ggplot(dat.merge.k9m3.merged, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.clsts2)
print(m.clsts2 + facet_wrap(~cluster))



# Show final before clean up  ---------------------------------------------

dat.bothmarks.pretty <- rbind(subset(dat.merge.k36.pretty, select = c(cell, umap1, umap2, umap1.scale, umap2.scale, umap1.shift, type, cluster, mark)),
                              subset(dat.merge.k9m3.merged, select = c(cell, umap1, umap2, umap1.scale, umap2.scale, umap1.shift, type, cluster, mark)))

dat.bothmarks.pretty$stage <- sapply(dat.bothmarks.pretty$cell, function(x) strsplit(x, "-")[[1]][[1]])
dat.bothmarks.pretty$stage <- factor(dat.bothmarks.pretty$stage, levels = c("E9p5", "E10", "E10p5", "E11p5"))

m.bothmarks.pretty <- ggplot(dat.bothmarks.pretty, aes(x = umap1.shift, y = umap2.scale, group = cell, color = cluster)) +
  geom_point(size = jsize) +
  ggtitle("K36 (left) and K9m3 (right) gastrulation E9.5, E10.5, E11.5", "color by celltype") +
  geom_path(alpha = 0.05) +
  scale_color_manual(values = cbPalette) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.bothmarks.pretty.stage <- ggplot(dat.bothmarks.pretty, aes(x = umap1.shift, y = umap2.scale, group = cell, color = stage)) +
  geom_point(size = jsize) +
  ggtitle("K36 (left) and K9m3 (right) gastrulation E9.5, E10.5, E11.5", "color by stage") +
  geom_path(alpha = 0.05) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



jcheck <- dat.bothmarks.pretty %>% filter(mark == "K9m3")

ggplot(jcheck, aes(x = umap1.shift, y = 1 * umap2.scale, color = cluster)) +
  geom_point() +
  facet_wrap(~cluster) +
  ggtitle("Double + single cells") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jcheck %>% mutate(cluster = ifelse(startsWith(x = cluster, prefix = "NeuralTubeNeuralProgs"), "NeuralTubeNeuralProgs", cluster)), aes(x = umap1.shift, y = 1 * umap2.scale, color = cluster)) +
  geom_point() +
  facet_wrap(~type) +
  ggtitle("Dbl and singles") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.bothmarks.pretty %>% filter(cluster != "NeuralTubeNeuralProgs2"), aes(x = umap1.shift, y = 1 * umap2.scale, group = cell, color = cluster)) +
  geom_point() +
  ggtitle("Double + single cells") +
  geom_path(alpha = 0.05) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.bothmarks.pretty %>% filter(cluster != "NeuralTubeNeuralProgs2" & cluster != "NeuralTubeNeuralProgs3"), aes(x = umap1.shift, y = 1 * umap2.scale, group = cell, color = cluster)) +
  geom_point() +
  ggtitle("Double + single cells") +
  geom_path(alpha = 0.05) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check on non-projected version  -----------------------------------------

# TODO

# Load the fits to check the bad cells ------------------------------------

inf.fits <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_outputs_objs/", jname, "/unmix_scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE.RData")
load(inf.fits, v=T)

fits.out <- act.repress.coord.lst
w.lst <- sapply(fits.out, function(x) x$w)

cell.vec <- unique(dat.merge.k9m3.dbl$cell)
names(cell.vec) <- cell.vec

# if louvains are now from clusters need eto rethink jcoord
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  jmax.LL <- jfit$ll.mat[jcoord]

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
  out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, logL = jmax.LL, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()

coords.dbl <- coords.dbl %>%
  rowwise() %>%
  mutate(celltype = clst.annots[[louv.act]])

m.grid <- ggplot(coords.dbl, aes(x = celltype, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")

print(m.grid)


dat.merge.k36.pretty.w <- left_join(dat.merge.k36.pretty %>% filter(type == "dbl"), coords.dbl)
dat.merge.k9m3.merged.w <- left_join(dat.merge.k9m3.merged %>% filter(type == "dbl"), coords.dbl)

ggplot(dat.merge.k36.pretty.w, aes(x = umap1, y = umap2, color = w)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge.k9m3.merged.w, aes(x = umap1, y = umap2, color = w)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.bothmarks.pretty.w <- subset(dat.bothmarks.pretty, type == "dbl") %>%
  left_join(., coords.dbl)

m.w <- ggplot(dat.bothmarks.pretty.w, aes(x = umap1.shift, y = umap2.scale, color = w, group = cell)) +
  geom_point(size = jsize) +
  ggtitle("K36 (left) and K9m3 (right) gastrulation E9.5, E10.5, E11.5", "w=FracReadsBelongingToK36") +
  geom_path(alpha = 0.05) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.density <- ggplot(dat.bothmarks.pretty.w, aes(x = w, fill = cluster)) +
  geom_density() +
  ggtitle("K36 (left) and K9m3 (right) gastrulation E9.5, E10.5, E11.5", "w=FracReadsBelongingToK36") +
  facet_wrap(~cluster) +
  scale_fill_manual(values = cbPalette) +
  xlab("Fraction of Reads Belonging to K36") +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")

# Write outputs -----------------------------------------------------------

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_demux_cleaned_", jname)
dir.create(outdir)
outf <- file.path(outdir, paste0("demux_cleaned_filtered_", jname, ".", Sys.Date(), ".filt2.spread_", jspread, ".dbl_only.txt"))
outf2 <- file.path(outdir, paste0("demux_cleaned_filtered_", jname, ".", Sys.Date(), ".filt2.spread_", jspread, ".single_and_dbl.txt"))
outpdf <- file.path(outdir, paste0("demux_cleaned_filtered_", jname, ".", Sys.Date(), ".filt2.spread_", jspread, ".pdf"))
fwrite(dat.bothmarks.pretty.w, file = outf, sep = "\t", quote = FALSE)
fwrite(dat.bothmarks.pretty, file = outf2, sep = "\t", quote = FALSE)

pdf(outpdf, useDingbats = FALSE)

print(m.bothmarks.pretty)
print(m.bothmarks.pretty.stage)
print(m.w)
print(m.density)

dev.off()


# Check k36me3 nonproj ----------------------------------------------------

# TODO

print(m.bothmarks.pretty)
print(m.bothmarks.pretty.stage)
print(m.w)
print(m.density)


