# Jake Yeung
# Date of Creation: 2021-07-11
# File: ~/projects/scChIX/analysis_scripts/4-analyze_scChIX_outputs_projections_cleaned_up2_from_genebodies_bins.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
# jsettings$spread <- 8


# Load matrix  ------------------------------------------------------------

jmark1 <- "K36"; jmark2 <- "K9m3"; jmarks <- c(jmark1, jmark2)
names(jmarks) <- jmarks

# jprefix <- "from_LDA_cleaned"
jprefix <- "from_cleaned_K36genebodies_K9m3bins_merged"
jsuffix <- "dbl_cleaned"
# jsuffix <- "dbl_k36_k9m3_cleaned"
infs <- lapply(jmarks, function(jmark){
    inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/", jprefix, "/", jsuffix, "/Gastru_Unmixed_DblMark.", jsuffix, ".", jmark, ".RData")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


# load mats

metamain <- file.path(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_", jprefix), jsuffix)
assertthat::assert_that(dir.exists(metamain))

dat.meta.lst <- lapply(jmarks, function(jmark){
  infmeta <- file.path(metamain, paste0("celltyping_output_filt.", jmark, ".2021-07-11.rds"))
  dat.meta <- readRDS(infmeta)
  return(dat.meta)
})


# get umap from original


# out.lda.predict, count.mat.proj, out.objs

dat.umap.merge.lst <- lapply(jmarks, function(jmarktmp){

  # jmarktmp <- jmarks[[1]]
  load(infs[[jmarktmp]], v=T)

  tm.orig <- posterior(out.objs$out.lda)
  umap.out <- umap(tm.orig$topics, config = jsettings)
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
  return(dat.umap.merge)
})

# ggplot(dat.umap.orig.annot, aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   theme_bw() +
#   ggtitle(jmarktmp) +
#   scale_color_manual(values = cbPalette) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.lst <- lapply(jmarks, function(jmarktmp){
  dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    facet_wrap(~type) +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

multiplot(m.lst[[1]], m.lst[[2]], cols = 1)

m.lst2 <- lapply(jmarks, function(jmarktmp){
  dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    # facet_wrap(~type) +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
multiplot(m.lst2[[1]], m.lst2[[2]], cols = 1)



# Split things up  --------------------------------------------------------


dat.merge.rbind <- bind_rows(dat.umap.merge.lst) %>%
  group_by(mark) %>%
  mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
         umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
         umap1.shift = ifelse(mark == "K36", umap1.scale - 5, umap1.scale + 5))

ggplot(dat.merge.rbind, aes(x = umap1.shift, y = -1 * umap2.scale, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.15) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge.rbind %>% filter(type == "dbl"), aes(x = umap1.shift, y = umap2.scale, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.15) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(dat.merge.rbind %>% filter(! (umap2.scale > -1 & umap1.shift < 0) ), aes(x = umap1.shift, y = umap2.scale, group = cell)) +
#   geom_point() +
#   geom_path(alpha = 0.05) +
#   theme_bw() +
#   theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Plot 2D -----------------------------------------------------------------


# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/double_stain_outputs"
# outrds <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".act_repress_coord_merged_lst.rds"))

hubprefix <- "/home/jyeung/hub_oudenaarden"
infrdata <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline/", jprefix, "/scchix_outputs.", jsuffix, "/unmix_scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE.RData"))
assertthat::assert_that(file.exists(infrdata))
load(infrdata, v=T)


# load(infrdata, v=T)

fits.out <- act.repress.coord.lst
w.lst <- sapply(fits.out, function(x) x$w)

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


library(ggforce)
m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

head(coords.dbl)

coords.dbl <- coords.dbl %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]],
         stage = gsub("E9", "E09", stage))

m.grid.bystage <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  facet_wrap(~stage) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.bystage)

m.grid.bystage <- ggplot(coords.dbl, aes(x = "a", y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  facet_wrap(~stage) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.bystage)

ggplot(coords.dbl, aes(x = stage, y = louv.repress)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")

coords.dbl.sum <- coords.dbl %>%
  group_by(louv.repress, stage) %>%
  summarise(ncount = length(louv.repress)) %>%
  ungroup() %>%
  mutate(nfrac = ncount / sum(ncount))

coords.dbl.sum.act <- coords.dbl %>%
  group_by(louv.act, stage) %>%
  summarise(ncount = length(louv.act)) %>%
  ungroup() %>%
  mutate(nfrac = ncount / sum(ncount))

ggplot(coords.dbl.sum, aes(x = stage, y = nfrac, color = louv.repress, group = louv.repress)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(coords.dbl.sum.act, aes(x = stage, y = nfrac, color = louv.act, group = louv.act)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Plot similarities onto 2D grid  -----------------------------------------

dat.merge.umap.dist.summary <- dat.merge.rbind %>%
  filter(type == "dbl") %>%
  dplyr::select(cell, umap1, umap2, mark) %>%
  reshape2::melt(., id.vars = c("cell", "mark"), measure.vars = c("umap1", "umap2")) %>%
  mutate(newvar = interaction(mark, variable, sep = "_")) %>%
  reshape2::dcast(data = ., formula = cell ~ newvar, value.var = "value")

coords.dbl.annot <- coords.dbl %>%
  left_join(., dat.merge.umap.dist.summary) %>%
  mutate(K36_avg = sqrt(K36_umap1 ^ 2 + K36_umap2 ^ 2),
         K9m3_avg = sqrt(K9m3_umap1 ^ 2 + K9m3_umap2 ^ 2))

actclsts <- unique(coords.dbl.annot$louv.act)
# jmarktest <- "K9m3"
jmarktest <- "K36"
jclst <- "cluster7"
for (jclst in actclsts){
  print(jclst)
  m.grid2 <- ggplot(coords.dbl.annot %>% filter(louv.act == jclst),
                    aes_string(x = "louv.act", y = "louv.repress", color = paste(jmarktest, "avg", sep = "_"))) +
    geom_point(alpha = 1, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
    theme_bw() +
    scale_color_viridis_c() +
    theme(aspect.ratio=2) +
    ggtitle(paste("Active cluster", jclst))
  print(m.grid2)
}

repclsts <- unique(coords.dbl.annot$louv.repress)
jmarktest <- "K9m3"
# jmarktest <- "K36"
# jclst <- "cluster7"
for (jclst in repclsts){
  print(jclst)
  m.grid2 <- ggplot(coords.dbl.annot %>% filter(louv.repress == jclst),
                    aes_string(x = "louv.act", y = "louv.repress", color = paste(jmarktest, "avg", sep = "_"))) +
    geom_point(alpha = 1, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
    theme_bw() +
    scale_color_viridis_c() +
    theme(aspect.ratio=0.5) +
    ggtitle(paste("Repress cluster", jclst))
  print(m.grid2)
}


# UMAP by stage -----------------------------------------------------------

dat.tmp <- dat.umap.merge.lst$K36  %>%
  rowwise() %>%
  mutate(stage = gsub("E9", "E09", strsplit(cell, split = "-")[[1]][[1]]))

ggplot(dat.tmp %>% filter(cluster != "na"), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_grid(cluster ~ stage) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.tmp %>% filter(cluster != "na"), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~stage) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.tmp %>% filter(cluster == "na"), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~stage) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


m.lst2 <- lapply(jmarks, function(jmarktmp){
  dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
  m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    # facet_wrap(~type) +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
multiplot(m.lst2[[1]], m.lst2[[2]], cols = 1)







