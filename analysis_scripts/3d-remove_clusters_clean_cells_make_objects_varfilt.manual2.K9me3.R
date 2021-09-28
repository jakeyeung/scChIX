# Jake Yeung
# Date of Creation: 2021-07-19
# File: ~/projects/scChIX/analysis_scripts/3d-remove_clusters_clean_cells_make_objects_varfilt.manual2.K9me3.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(scChIX)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123



# Load new LDA after cleaning  --------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jdate <- "2021-07-14"

jmarks <- c("K36", "K9m3", "K36-K9m3")
# jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

jprefix <- "var_filtered"
jdate <- "2021-07-16"
jquant <- "0.15"


# jcutoff.stringent <- exp(0.2)
# log.jcutoff.stringent <- log(jcutoff.stringent)

log.jcutoff.stringent.lst <- list(0.2, -0.5, -0.7)
names(log.jcutoff.stringent.lst) <- jmarks

bad.clsts <- list(c("cluster5"), c("cluster3"), c("cluster2"))
names(bad.clsts) <- jmarks

# Load LDA, meta, countmat objs  ----------------------------------------------------------

dname <- paste0(jprefix, "_", jquant, "_", jstr)

inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA", dname)

# inf.ldas <- lapply(jmarks, function(jmark){
#   inf.tmp <- file.path(inmain, paste0("lda_output_filt.", jmark, ".", jdate, ".rds"))
#   assertthat::assert_that(file.exists(inf.tmp))
#   return(inf.tmp)
# })


inf.countmats <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(inmain, paste0("countmat_output_filt.", jmark, ".", jdate, ".rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})


inf.metas <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(inmain, paste0("celltyping_output_filt.", jmark, ".", jdate, ".rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

# load datvarmat
jdate2 <- "2021-07-15"
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

inf.var <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2", paste0(jprefix, "_", jquant), jstr, paste0("intrachrom_var_outputs.", jmark, ".", jdate2, ".rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.var.lst <- lapply(inf.var, readRDS)


# print(dim(dat.var.lst[[jmark]]))
lapply(dat.var.lst, dim)


# jmark <- jmarks[[3]]
# m2 <- ggplot(dat.var.lst[[jmark]] %>% filter(cluster != "cluster1"), aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   theme_bw() +
#   ggtitle(jmark) +
#   scale_color_manual(values = cbPalette) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# print(m2)


m1.lst <- lapply(jmarks, function(jmark){
  m1 <- ggplot(dat.var.lst[[jmark]], aes(x = umap1, y = umap2, color = log(cell.var.within.sum.norm))) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_viridis_c(direction = -1) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m1)
})

m2.lst <- lapply(jmarks, function(jmark){
  m2 <- ggplot(dat.var.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m2)
})

m3.lst <- lapply(jmarks, function(jmark){
  m3 <- ggplot(dat.var.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    facet_wrap(~stage) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})





# Filter out bad clusters  ------------------------------------------------



# Filter dbl cells that are bad -------------------------------------------


jmarkdbl <- jmarks[[3]]
jname <- paste(jprefix, "manual", jstr, sep = "_")
infrdata <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_outputs_objs/", jname, "/unmix_scchix_inputs_clstr_by_celltype_", jmarkdbl, ".removeNA_FALSE.RData"))
assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)

fits.out <- act.repress.coord.lst
w.lst <- sapply(fits.out, function(x) x$w)

cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec

# if louvains are now from clusters need eto rethink jcoord
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


m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

# bad louvs from dbl
bad.pairs.act <- c("cluster1")

bad.cells.dbl <- subset(coords.dbl, louv.act %in% bad.pairs.act)$cell

m.grid.check <- ggplot(coords.dbl %>% filter(!cell %in% bad.cells.dbl), aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.check)



# Save outputs: countmat and meta ------------------------------------------------------------

countmats <- lapply(inf.countmats, readRDS)
datmetas <- lapply(inf.metas, readRDS)

outmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/var_filtered_manual2")
dir.create(outmain)
outdir <- file.path(outmain, jstr)
dir.create(outdir)

outpdf <- file.path(outdir, paste0("plots_var_filtering.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)
print(m1.lst)
print(m2.lst)
print(m3.lst)
m2.filt.lst <- lapply(jmarks, function(jmark){

  jcutoff.stringent <- exp(log.jcutoff.stringent.lst[[jmark]])

  m2.filt1 <- ggplot(dat.var.lst[[jmark]] %>% filter(cell.var.within.sum.norm >= jcutoff.stringent), aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  m2.filt2 <- ggplot(dat.var.lst[[jmark]] %>% mutate(is.good = cell.var.within.sum.norm >= jcutoff.stringent & !cell %in% bad.cells.dbl),
                     aes(x = umap1, y = umap2, color = is.good)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  m2.filt3 <- ggplot(dat.var.lst[[jmark]] %>% mutate(is.good = cell.var.within.sum.norm >= jcutoff.stringent & !cell %in% bad.cells.dbl), aes(x = log(cell.var.within.sum.norm))) +
    geom_density() +
    geom_vline(xintercept = log(jcutoff.stringent)) +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m2.filt1)
  print(m2.filt2)
  print(m2.filt3)
  return(m2.filt1)
})
m2.clstfilt.lst <- lapply(jmarks, function(jmark){
  m2.filt1 <- ggplot(dat.var.lst[[jmark]] %>%
                       filter(!cluster %in% bad.clsts[[jmark]]) %>%
                       filter(!cell %in% bad.cells.dbl),
                     aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})
print(m2.clstfilt.lst)
dev.off()

for (jmark in jmarks){

  jcutoff.stringent <- exp(log.jcutoff.stringent.lst[[jmark]])

  outcountmat <- file.path(outdir, paste0("countmat_var_filt.", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("metadata_var_filt.", jmark, ".", Sys.Date(), ".rds"))
  cells.keep <- subset(dat.var.lst[[jmark]], cell.var.within.sum.norm >= jcutoff.stringent & !cluster %in% bad.clsts[[jmark]] &  !cell %in% bad.cells.dbl)$cell
  cols.keep <- colnames(countmats[[jmark]]) %in% cells.keep
  countmat.tmp <- countmats[[jmark]][, cols.keep]
  print(dim(countmats[[jmark]]))
  print(dim(countmats[[jmark]][, cols.keep]))
  datmetas.tmp <- subset(datmetas[[jmark]], cell %in% colnames(countmat.tmp))
  saveRDS(countmat.tmp, file = outcountmat)
  saveRDS(datmetas.tmp, file = outmeta)
}

# check size
inf.check <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/var_filtered_manual/", jstr, "/countmat_var_filt.", jmarkdbl, ".2021-07-19.rds")
jcheck <- readRDS(inf.check)
print(dim(jcheck))

