# Jake Yeung
# Date of Creation: 2021-08-02
# File: ~/projects/scChIX/analysis_scripts/5-transfer_labels_from_scChIX.K9me3.check_K36.R
# description

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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load matrix  ------------------------------------------------------------

# jquants <- c("0.15", "0.2")
# jquants <- c("0.2")
# jquants <- c("0.15")
# jquants <- c("0.15")

# jdate <- "2021-07-16"
# jquant <- "0.15"

# jdate <- "2021-07-19"
# jquant <- "manual"

# jdate <- "2021-07-22"
# jquant <- "manual2noblood"

# jdate <- "2021-07-20"
# jquant <- "manual2"

jdate <- "2021-07-23"
jquant <- "manual2nocenter"

# jmark1 <- "K36"; jmark2 <- "K27"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")
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
# jsub <- subset(dat.merge.rbind, type == "dbl" & mark == jmarks[[2]] & umap2.scale < 2)
# jsub <- subset(dat.merge.rbind, type == "dbl" & mark == jmarks[[1]] & umap2.scale > -2)
# jsub <- subset(dat.merge.rbind, type == "dbl" & mark == jmarks[[2]])
# cell2lab <- hash::hash(jsub$cell, jsub$umap2.scale)
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



#
# # Transfer labels  --------------------------------------------------------
#
# # create KNN object
#
# # make into function?
# jmark <- jmarks[[1]]
#
# topics.mat1 <- tm.lst[[jmark]]$topics
# topics.mat2 <- out.lst.lst[[jmark]]$out.lda.predict$topics
#
# topics.mat.merge <- rbind(topics.mat1, topics.mat2)
#
# # run KNN
# umap.merged.out <- umap(topics.mat.merge, config = jsettings)
#
# # for dbl cells, get the nearest neighbor indices
# cells.all <- rownames(topics.mat.merge)
# names(cells.all) <- cells.all
#
# # get labels
# cell2cluster <- hash::hash(dat.merge.rbind$cell, dat.merge.rbind$cluster)
#
# dbl.cells <- rownames(out.lst.lst[[jmark]]$out.lda.predict$topics)
# names(dbl.cells) <- dbl.cells
#
# nn.indx.lst <- lapply(dbl.cells, function(jcell){
#   return(as.character(umap.merged.out$knn$indexes[jcell, ]))
# })
#
# # indx 2 cell
# indx2cell <- hash::hash(as.character(seq(cells.all)), cells.all)
#
# # get labels for each cell
# nn.cell.lst <- lapply(nn.indx.lst, function(indx.vec){
#   jvec <- sapply(indx.vec, function(x) AssignHash(x = x, jhash = indx2cell), USE.NAMES = FALSE)
#   names(jvec) <- jvec
#   return(jvec)
# })
#
# nn.clstr.lst <- lapply(nn.cell.lst, function(cell.vec){
#   sapply(cell.vec, function(x) AssignHash(x = x, jhash = cell2cluster))
# })
#
# dat.clstr.lst <- lapply(nn.clstr.lst, function(clstr.vec){
#   dat.clstr <- data.frame(cell = names(clstr.vec), cluster = clstr.vec, stringsAsFactors = FALSE)
# })
#
#
# dat.clstr.sum <- lapply(dbl.cells, function(jcell){
#   jdat <- dat.clstr.lst[[jcell]]
#   jdat.sum <- jdat %>%
#     group_by(cluster) %>%
#     summarise(ncells = length(cell),
#               nfrac = ncells / nrow(jdat)) %>%
#     filter(cluster != "na") %>%
#     arrange(desc(ncells)) %>%
#     mutate(cell = jcell)
#   if (nrow(jdat.sum) == 0){
#     data.frame(cluster = "na", ncells = nrow(jdat), nfrac = 1, cell = jcell, stringsAsFactors = FALSE)
#   } else {
#     return(jdat.sum[1, ])
#   }
# }) %>%
#   bind_rows()



# Do transfer but one cell at a time  -------------------------------------

CellToNearestNeighbors <- function(xvec, xmat, jmeth = "euclidean"){
  assertthat::assert_that(jmeth == "euclidean")
  xdist <- apply(xmat, 1, function(jrow) sum((xvec - jrow) ^ 2))
  xdist.sort <- sort(xdist, decreasing=FALSE)
  return(xdist.sort)
}

nnearest <- 30
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
  facet_wrap(~type) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Annotate clusters  ------------------------------------------------------



print(m.clsts)

# topic10 = cluster6 = "erythroid"
clst.annots <- list()
clst.annots["cluster6"] <- "Erythroid"
clst.annots["cluster1"] <- "Endothelial"
clst.annots["cluster3"] <- "Epithelial"
clst.annots["cluster2"] <- "Neurons"
clst.annots["cluster11"] <- "ConnectiveTissueProg"
clst.annots["cluster9"] <- "WhiteBloodCells"
clst.annots["cluster4"] <- "NeuralTubeNeuralProgs"
clst.annots["cluster10"] <- "Stromal"
clst.annots["cluster7"] <- "SchwannCellPrecusor"

# merge cluster5 and cluster8 into neural progs
clst.annots["cluster5"] <- "NeuralTubeNeuralProgs2"
clst.annots["cluster8"] <- "NeuralTubeNeuralProgs3"

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

dat.merge.k36.prettyy <- dat.merge.k36 %>%
  dplyr::rename(cluster.louvain = cluster,
                cluster = celltype)

# Write meta data ---------------------------------------------------------

outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix"
outdir <- file.path(outmain, jname)
dir.create(outdir)

outf <- file.path(outdir, paste0("celltyping_", jmark, "_first_try.", Sys.Date(), ".txt"))
outpdf <- file.path(outdir, paste0("celltyping_", jmark, "_first_try.", Sys.Date(), ".pdf"))
fwrite(x = dat.merge.k36, file = outf, sep = "\t")

pdf(outpdf, useDingbats = FALSE)
  print(m.clsts)
  print(m.annot)
dev.off()
