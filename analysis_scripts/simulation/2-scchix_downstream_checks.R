# Jake Yeung
# Date of Creation: 2021-09-19
# File: ~/projects/scChIX/analysis_scripts/simulation/2-scchix_downstream_checks.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(dplyr)

library(topicmodels)

hubprefix <- "/home/jyeung/hub_oudenaarden"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load sim params ---------------------------------------------------------

inf.sim <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_inputs/countmats/ATAC_simulator_params_and_outputs.RData")
load(inf.sim, v=T)


# Load unmixing mat -------------------------------------------------------

inf.prob <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_outputs/scchix_unmixing_downstream/scchix_inputs_clstr_by_celltype-prob_mat.mark1-mark2_to_mark1.txt")
dat.prob <- fread(inf.prob)
mat.prob <- as.matrix(subset(dat.prob, select = -V1))
rownames(mat.prob) <- dat.prob$V1

plot(density(mat.prob[, 1]))



# Load unmixing  ----------------------------------------------------------



# Check outputs -----------------------------------------------------------

inf.input <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_outputs/scchix_inputs_objs/scchix_inputs_clstr_by_celltype_mark1-mark2.removeNA_FALSE.RData")
load(inf.input, v=T)

inf.output <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_mark1-mark2.RData")
load(inf.output, v=T)

fits.out <- act.repress.coord.lst
w.lst <- sapply(fits.out, function(x){
  return(x$w)
})
plot(hist(w.lst))  # more active than repressive? Why is that?

plot(density(w.lst))

out.dat <- data.frame(cell = names(w.lst), w = w.lst, stringsAsFactors = FALSE)

# bimodal???
ggplot(out.dat, aes(x = w)) + geom_histogram() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# get active, repress index for each cell
cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  mark1.clstname <- rownames(jfit$ll.mat)[[jcoord[[1]]]]
  mark2.clstname <- colnames(jfit$ll.mat)[[jcoord[[2]]]]
  out.dat <- data.frame(cell = jcell, louv.act = mark1.clstname, louv.repress = mark2.clstname, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()

library(ggforce)
ggplot(coords.dbl, aes(x = as.character(louv.act), y = as.character(louv.repress))) +
  geom_point(position = position_jitternormal(sd_x = 0.1, sd_y = 0.1), alpha = 0.25) +
  theme_bw() +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Active Clusters") + ylab("Repressed Clusters") +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")


# Plot UMAP  --------------------------------------------------------------

inf.meta.dbl <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_outputs/objs_from_LDA/celltyping_output_filt.mark1-mark2.rds")
dat.meta.dbl <- readRDS(inf.meta.dbl) %>%
  rowwise() %>%
  mutate(ctype = strsplit(cell, "_")[[1]][[3]])


ggplot(dat.meta.dbl, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~ctype) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

coords.dbl.annot <- left_join(coords.dbl, subset(dat.meta.dbl, select = c(cell, umap1, umap2, ctype)))

ggplot(coords.dbl.annot, aes(x = as.character(louv.act), y = as.character(louv.repress), color = ctype)) +
  geom_point(position = position_jitternormal(sd_x = 0.1, sd_y = 0.1), alpha = 0.25) +
  theme_bw() +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Active Clusters") + ylab("Repressed Clusters") +
  ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")



# Check other  ------------------------------------------------------------

jmarks <- c("mark1", "mark2", "mark1-mark2"); names(jmarks) <- jmarks

inf.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/simulation_data/snakemake_outputs/objs_from_LDA/celltyping_output_filt.", jmark, ".rds"))
  return(inf.meta)
})

dat.meta.lst <- lapply(jmarks, function(jmark){
  jinf <- inf.meta.lst[[jmark]]
  dat.meta.tmp <- readRDS(jinf) %>%
  rowwise() %>%
  mutate(ctype = substr(cell, start = nchar(cell), nchar(cell)),
         mark = jmark)
})

m.meta.ctype <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

m.meta.cluster <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

JFuncs::multiplot(m.meta.ctype[[1]], m.meta.ctype[[2]], m.meta.ctype[[3]], cols = 3)
JFuncs::multiplot(m.meta.cluster[[1]], m.meta.cluster[[2]], m.meta.cluster[[3]], cols = 3)



# Connect the UMAPs  ------------------------------------------------------


inf1 <- file.path(hubprefix, "/jyeung/data/dblchic/simulation_data/snakemake_outputs/projection_output.mark1.RData")
inf2 <- file.path(hubprefix, "/jyeung/data/dblchic/simulation_data/snakemake_outputs/projection_output.mark2.RData")

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
jsettings1$spread <- 1

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

dat.umap.out1.merged <- rbind(dat.umap.out1, dat.umap.out1.pred) %>%
  mutate(mark = "mark1")
dat.umap.out2.merged <- rbind(dat.umap.out2, dat.umap.out2.pred) %>%
  mutate(mark = "mark2")

ggplot(dat.umap.out1.merged, aes(x = umap1, y = umap2)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.out2.merged, aes(x = umap1, y = umap2)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Join umaps --------------------------------------------------------------

dat.umap.merged <- rbind(dat.umap.out1.merged, dat.umap.out2.merged)

dat.final.annots <- dat.meta.lst %>%
  bind_rows()

dat.final.annots.sub <- subset(dat.final.annots, select = c(cell, cluster, ctype))

dat.umap.joined <- left_join(dat.umap.merged, dat.final.annots.sub) %>%
  group_by(mark) %>%
  mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
         umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
         umap1.scale.shift = ifelse(mark == "mark1", umap1.scale - 3, umap1.scale + 3))


ggplot(dat.umap.joined, aes(x = umap1, y = umap2, color = ctype)) +
  geom_point() +
  facet_wrap(~mark, nrow = 1) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = ctype, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.02) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = ctype, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.02) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  facet_wrap(~cluster) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")




# Verify we can infer mutual and overlapping bins  ------------------------

inf.mats <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/simulation_data/snakemake_inputs/countmats/countmat_var_filt.", jmark, ".rds"))
})
mats <- lapply(inf.mats, readRDS)


bin.dat.A <- ctype.sim.counts.lst$A$bin.data
bin.dat.A <- bin.dat.A[gtools::mixedorder(bin.dat.A$Bin), ]

mat.prob.A <- mat.prob[, grepl("_A$", colnames(mat.prob))]
mat.A.dbl <- mats$`mark1-mark2`[, grepl("_A$", colnames(mats$`mark1-mark2`))]

plot(density(colSums(mat.A.dbl)))

plot(density(mat.prob.A[1, ]))
plot(density(mat.prob.A[2, ]))

plot(density(rowMeans(mat.prob.A)))
plot(hist(rowMeans(mat.prob.A), breaks = 100))


# compare with ground truth
jnbins <- nrow(bin.dat.A)
frac.swap <- ctype.params.lst$A$frac.mutual.excl / 2
bottom.bins <- (bin.dat.A %>% dplyr::arrange(BinMean))$Bin[1:(jnbins * frac.swap)]
top.bins <- (bin.dat.A %>% dplyr::arrange(desc(BinMean)))$Bin[1:(jnbins * frac.swap)]


bin.means.A1 <- data.frame(mark = "mark1", BinMean = bin.dat.A$BinMean, stringsAsFactors = FALSE)
rownames(bin.means.A1) <- bin.dat.A$Bin

bin.means.A2 <- data.frame(mark = "mark2", BinMean = bin.dat.A$BinMean, stringsAsFactors = FALSE)
rownames(bin.means.A2) <- bin.dat.A$Bin

SwapBins <- function(ctype1.sim.counts.b, top.bins, bottom.bins, frac.swap){
  ctype1.sim.counts.b.swapped <- ctype1.sim.counts.b
  ctype1.sim.counts.b.swapped[bottom.bins, ] <- ctype1.sim.counts.b[top.bins, ]
  ctype1.sim.counts.b.swapped[top.bins, ] <- ctype1.sim.counts.b[bottom.bins, ]
  return(ctype1.sim.counts.b.swapped)
}

bin.means.A2 <- SwapBins(ctype1.sim.counts.b = bin.means.A2, top.bins = top.bins, bottom.bins = bottom.bins, frac.swap = frac.swap)

bin.means.A1$Bin <- rownames(bin.means.A1)
bin.means.A2$Bin <- rownames(bin.means.A2)

bin.means.A.merged <- left_join(bin.means.A1, bin.means.A2, by = "Bin") %>%
  rowwise() %>%
  mutate(BinMean.dbl = (BinMean.x) / (BinMean.x + BinMean.y))

plot(hist(bin.means.A.merged$BinMean.dbl, breaks = 100))


plot(bin.means.A.merged$BinMean.dbl, rowMeans(mat.prob.A))




