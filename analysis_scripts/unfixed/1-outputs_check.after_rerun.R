# Jake Yeung
# Date of Creation: 2021-09-15
# File: ~/projects/scChIX/analysis_scripts/unfixed/1-outputs_check.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(ggforce)



# Load projects -----------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf1 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/projection_output.K4m1.RData")
inf2 <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/projection_output.K27m3.RData")

load(inf1, v=T)
out.lda.predict1 <- out.lda.predict
count.mat.proj1 <- count.mat.proj
out.objs1 <- out.objs
load(inf2, v=T)
out.lda.predict2 <- out.lda.predict
count.mat.proj2 <- count.mat.proj
out.objs2 <- out.objs


# Plot outputs ------------------------------------------------------------

jsettings1 <- umap.defaults
jsettings1$n_neighbors <- 30
jsettings1$min_dist <- 0.1
jsettings1$random_state <- 123
jsettings1$spread <- 7

jsettings2 <- umap.defaults
jsettings2$n_neighbors <- 30
jsettings2$min_dist <- 0.1
jsettings2$random_state <- 123
jsettings2$spread <- 1

umap.out1 <- umap(posterior(out.objs1$out.lda)$topics, config = jsettings1)
umap.out2 <- umap(posterior(out.objs2$out.lda)$topics, config = jsettings2)


umap.out1.pred <- predict(umap.out1, out.lda.predict1$topics)
umap.out2.pred <- predict(umap.out2, out.lda.predict2$topics)

dat.umap.out1.pred <- data.frame(cell = rownames(umap.out1.pred), umap1 = umap.out1.pred[, 1], umap2 = umap.out1.pred[, 2], experi = "double", stringsAsFactors = FALSE)
dat.umap.out2.pred <- data.frame(cell = rownames(umap.out2.pred), umap1 = umap.out2.pred[, 1], umap2 = umap.out2.pred[, 2], experi = "double", stringsAsFactors = FALSE)

dat.umap.out1 <- data.frame(cell = rownames(umap.out1$layout), umap1 = umap.out1$layout[, 1], umap2 = umap.out1$layout[, 2], experi = "single")
dat.umap.out2 <- data.frame(cell = rownames(umap.out2$layout), umap1 = umap.out2$layout[, 1], umap2 = umap.out2$layout[, 2], experi = "single")

# umap.out1.merged <- rbind(dat.umap.out1, umap
# umap.out2.merged <- rbind(dat.umap.out2, umap.out2.pred)
# dat.umap.out1.merged <- data.frame(cell = rownames(umap.out1.merged), umap1 = umap.out1.merged[, 1], umap2 = umap.out1.merged[, 2], stringsAsFactors = FALSE) %>%
#   mutate(mark = "K4m1")
# dat.umap.out2.merged <- data.frame(cell = rownames(umap.out2.merged), umap1 = umap.out2.merged[, 1], umap2 = umap.out2.merged[, 2], stringsAsFactors = FALSE) %>%
#   mutate(mark = "K27m3")

dat.umap.out1.merged <- rbind(dat.umap.out1, dat.umap.out1.pred) %>%
  mutate(mark = "K4m1")
dat.umap.out2.merged <- rbind(dat.umap.out2, dat.umap.out2.pred) %>%
  mutate(mark = "K27m3")

ggplot(dat.umap.out1.merged, aes(x = umap1, y = umap2)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.out2.merged, aes(x = umap1, y = umap2)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load metas --------------------------------------------------------------

inf.rdata <- file.path(hubprefix, "jyeung/data/dblchic/from_rstudio/primetime/unfixed_louvain2/BM_UnfixedLouvain2.FinalCellClusterTable.2020-03-21.RData")
load(inf.rdata, v=T)

dat.umap.merged <- rbind(dat.umap.out1.merged, dat.umap.out2.merged)

dat.final.annots.sub <- subset(dat.final.annots, select = c(cell, experi, mark, cluster))

dat.umap.joined <- left_join(dat.umap.merged, dat.final.annots.sub) %>%
  group_by(mark) %>%
  mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
         umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
         umap1.scale.shift = ifelse(mark == "K4m1", umap1.scale - 3, umap1.scale + 3))


ggplot(dat.umap.joined, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~mark, nrow = 1) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = cluster, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.02) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = cluster, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.02) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  facet_wrap(~cluster) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

#
# # Check outputs -----------------------------------------------------------
#
# inf.input <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/scchix_inputs_objs/scchix_inputs_clstr_by_celltype_K4m1-K27m3.removeNA_FALSE.RData")
# load(inf.input, v=T)
#
# inf.output <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_K4m1-K27m3.RData")
# load(inf.output, v=T)
#
# fits.out <- act.repress.coord.lst
# w.lst <- sapply(fits.out, function(x){
#   return(x$w)
# })
# plot(hist(w.lst))  # more active than repressive? Why is that?
#
# out.dat <- data.frame(cell = names(w.lst), w = w.lst, stringsAsFactors = FALSE)
# # out.dat$experi <- sapply(as.character(out.dat$cell), function(x) paste(strsplit(x, "_")[[1]][1:2], collapse = "_"))
#
# dat.umap.dbl.merge <- left_join(dat.umap.joined, out.dat)
#
# m.w <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y = umap2, color = w)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_viridis_c()
#
# print(m.w)
#
# # bimodal???
# ggplot(out.dat, aes(x = w)) + geom_histogram() +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
#
# # get active, repress index for each cell
# cell.vec <- names(fits.out)
# names(cell.vec) <- cell.vec
# coords.dbl <- lapply(cell.vec, function(jcell){
#   jfit <- fits.out[[jcell]]
#   jweight <- fits.out[[jcell]]$w
#   p.mat <- SoftMax(jfit$ll.mat)
#   jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
#   jmax <- max(p.mat)
#   mark1.clstname <- rownames(jfit$ll.mat)[[jcoord[[1]]]]
#   mark2.clstname <- colnames(jfit$ll.mat)[[jcoord[[2]]]]
#   out.dat <- data.frame(cell = jcell, louv.act = mark1.clstname, louv.repress = mark2.clstname, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
#   return(out.dat)
# }) %>%
#   bind_rows()
#
#
# ggplot(coords.dbl, aes(x = louv.act, y = louv.repress)) +
#   geom_point(alpha = 0.02) +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
#
#
# logpcutoff <- log(0.9999)
# coords.dbl.stringent <- subset(coords.dbl, lnprob > logpcutoff)
#
# # count all pairs
# dbl.pairs <- paste(coords.dbl.stringent$louv.act, coords.dbl.stringent$louv.repress, sep = "_")
# dbl.pairs.counts <- sort(table(dbl.pairs), decreasing = TRUE)
# dbl.pairs.counts.filt <- dbl.pairs.counts[which(dbl.pairs.counts >= 10)]
#
# plot(hist(dbl.pairs.counts, breaks = 100))
#
#
# coords.dbl.annots <- left_join(coords.dbl, dat.umap.joined)
#
# ggplot(coords.dbl.annots, aes(x = as.character(louv.act), y = as.character(louv.repress), color = cluster)) +
#   geom_point(position = position_jitternormal(sd_x = 0.1, sd_y = 0.1), alpha = 0.25) +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   xlab("Active Clusters") + ylab("Repressed Clusters") +
#   ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")
#
#
# # what's the purity for each pair?
# jsum <- coords.dbl.annots %>%
#   mutate(actxrep = interaction(louv.act, louv.repress)) %>%
#   group_by(actxrep, cluster) %>%
#   summarise(ncells = length(cell)) %>%
#   group_by(actxrep) %>%
#   mutate(nfrac = ncells / sum(ncells))
#
# jentrop <- jsum %>%
#   group_by(actxrep) %>%
#   summarise(H = -sum(nfrac * log(nfrac)),
#             ncells = sum(ncells)) %>%
#   arrange(desc(H))
#
#
# dat.final.annots.k27 <- subset(dat.final.annots, mark == "K27m3" & experi == "single") %>%
#   dplyr::rename(celltype = cluster) %>%
#   dplyr::select(cell, celltype)
#
# dat.final.annots.k4m1 <- subset(dat.final.annots, mark == "K4m1" & experi == "single") %>%
#   dplyr::rename(celltype = cluster) %>%
#   dplyr::select(cell, celltype)
#
#
# # plot cluster names onto umap
#
# m1.k4 <- ggplot(dat.louv.lst$K4m1 %>% left_join(., dat.final.annots.k4m1), mapping = aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# m2.k4 <- ggplot(dat.louv.lst$K4m1 %>% left_join(., dat.final.annots.k4m1), mapping = aes(x = umap1, y = umap2, color = celltype)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
# JFuncs::multiplot(m1.k4, m2.k4, cols = 2)
#
#
# m1 <- ggplot(dat.louv.lst$K27m3 %>% left_join(., dat.final.annots.k27), mapping = aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# m2 <- ggplot(dat.louv.lst$K27m3 %>% left_join(., dat.final.annots.k27), mapping = aes(x = umap1, y = umap2, color = celltype)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
# JFuncs::multiplot(m1, m2, cols = 2)
#
# ggplot(coords.dbl.annots, aes(x = umap1, y = umap2)) +
#   geom_point() +
#   facet_wrap(~mark) +
#   theme_bw() +
#   scale_color_manual(values = cbPalette, na.value = "grey85") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#
#
# # Filter out cells  -------------------------------------------------------
#
# clst.remove <- c("cluster6")
#
# cells.remove.single <- subset(dat.louv.lst$K4m1, cluster == "cluster6")$cell
# cells.remove.dbl <- subset(coords.dbl.annots, louv.act == "cluster6")$cell
#
# cells.remove <- c(cells.remove.single, cells.remove.dbl)
#
# cells.keep <- dat.final.annots$cell[!dat.final.annots$cell %in% cells.remove]
#
# # Load mats again ---------------------------------------------------------
#
# indir.mats <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_inputs/countmats.first_try")
# inf.mats <- list.files(indir.mats, pattern = "*.rds", full.names = TRUE)
# names(inf.mats) <- basename(unlist(inf.mats))
# jnames <- names(inf.mats)
#
# mats.input <- lapply(inf.mats, function(jinf){
#   readRDS(jinf)
# })
#
# mats.filt <- lapply(mats.input, function(jmat){
#   cols.keep <- colnames(jmat) %in% cells.keep
#   jmat.out <- jmat[, cols.keep]
# })
#
# # Write outputs -----------------------------------------------------------
#
# outdir <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/snakemake_10kb_pipeline/snakemake_inputs/countmats")
# for (jname in jnames){
#   outf <- file.path(outdir, jname)
#   saveRDS(mats.filt[[jname]], file = outf)
# }


