# Jake Yeung
# Date of Creation: 2021-04-15
# File: ~/projects/scChIX/interactive_scripts/1-subset_cells_for_fitting.R
#



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(hash)
library(umap)
library(ggforce)

# Constants ---------------------------------------------------------------

MakeNamesCells <- FALSE  # in this analysis the dots have returned to dashes. So no need to use make.names

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jctypes.order <- c("Granulocytes", "Bcells", "NKcells")

jmarks <- c("K27m3", "K9m3")
jmarks.dbl <- c("K27m3xK9m3")
names(jmarks) <- jmarks

jmarks.all <- c(jmarks, jmarks.dbl)
names(jmarks.all) <- jmarks.all

jsuffix <- ".FewerTopics2"
jprefix.base <- "TopBins_autosomesOnly"
jprefix <- paste0(jprefix.base, jsuffix, "_KeepAllPlates_clstrbytopics")
jprefix.input <- jprefix
jprefix.output <- paste0(jprefix, "_0_1")
jprefix.ldaoutput <- paste0(jprefix, "_0_1_RemoveBadMixings")

# indir.clsts <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_downstream/clusterbytopics_TopBins_autosomesOnly", jsuffix)
# assertthat::assert_that(dir.exists(indir.clsts))

inf.input <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/double_staining_input/", jprefix.input, "/", jprefix, "_", jmarks.dbl, ".removeNA_TRUE.RData")
inf.output <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/double_staining_output/", jprefix.output, "/unmix_", jprefix, "_", jmarks.dbl, ".removeNA_TRUE.RData")
indir.lda <-  paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/projections_LDA_outputs/afterUnmixing_", jprefix.ldaoutput, "_K27m-K9m3")

assertthat::assert_that(file.exists(inf.input))
assertthat::assert_that(file.exists(inf.output))
assertthat::assert_that(dir.exists(indir.lda))

# manual ordering (see pdf of unmixing downstream for which louvains to match) ------

jsep <- "_"
louv.act.ordering.all <- c("topic16", "topic3", "topic9", "topic25", "topic11")
louv.repress.ordering.all <- c("topic12", "topic28", "topic1", "topic20")

jtop.pairs <- list(c(louv.act.ordering.all[[1]], louv.repress.ordering.all[[1]]),
                   c(louv.act.ordering.all[[2]], louv.repress.ordering.all[[2]]),
                   c(louv.act.ordering.all[[3]], louv.repress.ordering.all[[2]]),
                   c(louv.act.ordering.all[[4]], louv.repress.ordering.all[[3]]),
                   c(louv.act.ordering.all[[5]], louv.repress.ordering.all[[4]]))
jtop.pairs.names <- sapply(jtop.pairs, function(jtop.pair) paste(jtop.pair, collapse = jsep))
names(jtop.pairs) <- jtop.pairs.names

louv.act.ordering <- unique(louv.act.ordering.all)
louv.repress.ordering <- unique(louv.repress.ordering.all)

# Load annots --------------------------------------------------------------

if (MakeNamesCells){
  annots.dat <- LoadCellAnnotsEtOH.MorePlates()
} else {
  annots.dat <- LoadCellAnnotsEtOH.MorePlates() %>%
    dplyr::rename(cell.makenames = cell,
                  cell = cell.orig)
}

# Load data  --------------------------------------------------------------

load(inf.input, v=T)
load(inf.output, v=T)
# Make grid ---------------------------------------------------------------


fits.out <- act.repress.coord.lst

w.lst <- sapply(fits.out, function(x) x$w)

# remove 0.01 or 0.99
cells.remove.i <- which(w.lst >= 0.99 | w.lst <= 0.01)
if (length(cells.remove.i) > 0){
  cells.remove <- names(w.lst)[cells.remove.i]
  fits.out[[cells.remove]] <- NULL
}


# if louvains are now from clusters need eto rethink jcoord
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


coords.dbl.annots <- left_join(coords.dbl, annots.dat)

coords.dbl.annots$louv.act <- factor(coords.dbl.annots$louv.act, levels = louv.act.ordering)
coords.dbl.annots$louv.repress <- factor(coords.dbl.annots$louv.repress, levels = louv.repress.ordering)

m.grid <- ggplot(coords.dbl.annots, aes(x = louv.act, y = louv.repress, color = celltype)) +
  geom_point(position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")


print(m.grid)
head(coords.dbl.annots)


# get some granus
cells.g <- subset(coords.dbl.annots, louv.act == "topic16" & louv.repress == "topic12" & celltype == "Granulocytes")$cell
cells.b <- subset(coords.dbl.annots, louv.act %in% c("topic3", "topic9")  & louv.repress == "topic28" & celltype == "Bcells")$cell
cells.n <- subset(coords.dbl.annots, louv.act %in% c("topic25", "topic11")  & louv.repress %in% c("topic1", "topic20") & celltype == "NKcells")$cell

keepn <- 30
keepn.i <- keepn + 1
# load mat
load("/home/jyeung/projects/scChIX/interactive_scripts/ObjsForFitting.RData", v=T)
cuts.total <- sort(colSums(count.mat.doubleincub), decreasing = TRUE)

cuts.cells.g.subset <- cuts.total[cells.g][2:keepn.i]
cuts.cells.b.subset <- cuts.total[cells.b][2:keepn.i]
cuts.cells.n.subset <- cuts.total[cells.n][2:keepn.i]
cells.g.subset <- names(cuts.cells.g.subset)
cells.b.subset <- names(cuts.cells.b.subset)
cells.n.subset <- names(cuts.cells.n.subset)

jsubset <- cells.g.subset
m.grid.subset <- ggplot(coords.dbl.annots %>% filter(cell %in% jsubset), aes(x = louv.act, y = louv.repress, color = celltype)) +
  geom_point(position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.subset)

count.mat.subset <- count.mat.doubleincub[, ]
jsubset <- cells.b.subset

m.grid.subset <- ggplot(coords.dbl.annots %>% filter(cell %in% jsubset), aes(x = louv.act, y = louv.repress, color = celltype)) +
  geom_point(position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.subset)


jsubset <- cells.n.subset
m.grid.subset <- ggplot(coords.dbl.annots %>% filter(cell %in% jsubset), aes(x = louv.act, y = louv.repress, color = celltype)) +
  geom_point(position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
  print(m.grid.subset)


count.mat.subset <- count.mat.doubleincub[, c(cells.g.subset, cells.b.subset, cells.n.subset)]
load("/home/jyeung/projects/scChIX/interactive_scripts/TopBins_autosomesOnly.FewerTopics2_KeepAllPlates_clstrbytopics_K27m3xK9m3.removeNA_TRUE.RData", v=T)
rkeep <- rownames(count.dat$counts)
jsub <- count.mat.subset[rkeep, ]
jsub.allcells <- count.mat.doubleincub[rkeep, ]
count.mat.dbl.subset <- jsub
count.mat.dbl.allcells <- jsub.allcells

outf.subset <- "/home/jyeung/projects/scChIX/data/RawDblCountMatSubset.RData"
outf.allcells <- "/home/jyeung/projects/scChIX/data/RawDblCountMatAllCells.RData"
outf.models <- "/home/jyeung/projects/scChIX/data/ModelFromSingleForFitting.RData"
outf2 <- "/home/jyeung/projects/scChIX/data/GroundTruthLabels.RData"
save(count.mat.dbl.subset, file = outf.subset)
save(count.mat.dbl.allcells, file = outf.allcells)
save(dat.impute.repress.lst, dat.impute.active, file = outf.models)
save(annots.dat, file = outf2)



