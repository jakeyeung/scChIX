# Jake Yeung
# Date of Creation: 2021-10-22
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/8-infer_counts_across_celltypes.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scChIX)
library(scchicFuncs)

hubprefix <- "/Users/yeung/hub_oudenaarden"

outpdf <- paste0("/Users/yeung/data/dblchic/gastrulation/H3K36me3_H3K9me3_downstream_plots/k36_k9_ratio_calculation.", Sys.Date(), ".pdf")

# Load prob mat  ----------------------------------------------------------


# 1GB so copy to disk first
# inf.mat <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/scchix_unmixing_downstream/scchix_inputs_clstr_by_celltype-prob_mat.K36-K9m3_to_K36.txt")
inf.mat <- "/Users/yeung/Dropbox/from_cluster/dblchic/gastrulation_objs/K36_K9m3_K36-K9m3/scchix_inputs_clstr_by_celltype-prob_mat.K36-K9m3_to_K36.txt"
prob.mat <- fread(inf.mat)
rnames <- prob.mat$V1
prob.mat <- as.matrix(subset(prob.mat, select = -V1))
rownames(prob.mat) <- rnames

# Load meta  --------------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/dat_meta_ordered_colorcoded.rds")
dat.meta.lst <- readRDS(inf.meta)
dat.meta.long <- dat.meta.lst %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(plate = scchicFuncs::ClipLast(x = cell, jsep = "_"))

plates <- unique(sapply(dat.meta.long$cell, function(x) scchicFuncs::ClipLast(x = x, jsep = "_")))
# dbl plates
dbl.plates <- plates[grepl("K36.*K9", plates)]

dat.meta.filt <- subset(dat.meta.long, mark == "K36" & plate %in% dbl.plates)

prob.avg.bycell <- colMeans(prob.mat)
dat.prob.avg <- data.frame(cell = names(prob.avg.bycell), prob = prob.avg.bycell, stringsAsFactors = FALSE) %>%
  left_join(., dat.meta.filt)

inf.meta.cleaned <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/metadata_cleaned/metadata_cleaned.2021-10-04.rds")
dat.meta.cleaned <- readRDS(inf.meta.cleaned)

# pdf(file = file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/prob_summaries/dat_prob_avg_plots.", Sys.Date(), ".pdf")), useDingbats = FALSE)

# plot on the UMAP

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette.coefficients <- cbPalette[2:length(cbPalette)]

dat.umap.prob <- left_join(subset(dat.prob.avg, select = c(cell, prob)), dat.meta.cleaned) %>%
  filter(!is.na(prob))

ggplot(dat.umap.prob, aes(x = umap1, y = umap2, color = log2(prob / (1 - prob)), group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.05) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


ggplot(dat.prob.avg, aes(x = prob, fill = cluster)) + 
  geom_density() + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~cluster, nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.prob.avg, aes(x = prob, fill = colorcode)) + 
  geom_density(alpha = 0.25) + 
  scale_fill_identity() + 
  # facet_wrap(~cluster, nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# get avg for erythroid
eryth.avg <- median(subset(dat.prob.avg, cluster == "Erythroid")$prob)
ggplot(dat.prob.avg, aes(x = prob, fill = colorcode)) + 
  geom_density(alpha = 1) + 
  geom_vline(xintercept = eryth.avg, linetype = "dotted") + 
  scale_fill_identity() + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob)) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / (K36 + K9)") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = log2(prob / (1 - prob)))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("log2(K36 / K9)") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob / (1 - prob))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / K9") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# dev.off()


# Can we know if K36 or K9me3 is higher from single counts?  --------------

# load raw mats? 
jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks

infs.mat <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_inputs/countmats/countmat_var_filt.", jmark, ".rds"))
  return(inf)
})

mats.lst <- lapply(infs.mat, function(jinf){
  readRDS(jinf)
})

dat.cells.sum.lst <- lapply(jmarks, function(jmark){
  jmat <- mats.lst[[jmark]]
  data.frame(totalcuts = colSums(jmat), cell = colnames(jmat), stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.lst[[jmark]]) %>%
    rowwise() %>%
    mutate(plate = ClipLast(x = cell, jsep = "_"))
})

dat.cells.sum.byclstr.lst <- lapply(dat.cells.sum.lst, function(jdat){
  jdat %>%
    group_by(cluster) %>%
    summarise(mean.linear = mean(totalcuts), 
              mean.log2 = mean(log2(totalcuts)))
})

ggplot(dat.cells.sum.lst$K36, aes(x = umap1, y = umap2, color = log2(totalcuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.cells.sum.lst$K36, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K36, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  facet_wrap(~plate) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  facet_wrap(~plate) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = umap1, y = umap2, color = log2(totalcuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Combine some sparse plates then plot the fold changes -------------------

print(unique(dat.cells.sum.lst$K9m3$plate))

plates.to.combine <- list("E10p5-B6C-K9m3-190409-1" = "E10p5-B6C-K9m3-190409-1x2", 
                          "E10p5-B6C-K9m3-190409-2" = "E10p5-B6C-K9m3-190409-1x2", 
                          "E10p5-B6C-K9m3-190409-3" = "E10p5-B6C-K9m3-190409-3x4",
                          "E10p5-B6C-K9m3-190409-4" = "E10p5-B6C-K9m3-190409-3x4", 
                          "E11p5-CB6-K9m3-190618-1" = "E11p5-CB6-K9m3-190618-1x2", 
                          "E11p5-CB6-K9m3-190618-2" = "E11p5-CB6-K9m3-190618-1x2", 
                          "E11p5-CB6-K9m3-190618-3" = "E11p5-CB6-K9m3-190618-3x4", 
                          "E11p5-CB6-K9m3-190618-4" = "E11p5-CB6-K9m3-190618-3x4", 
                          "E9p5-CB6-K9m3-190409-1" = "E10p5-B6C-K9m3-190409-1x2", 
                          "E9p5-CB6-K9m3-190409-2" = "E10p5-B6C-K9m3-190409-1x2")

dat.cells.sum.lst$K9m3$platemerged <- sapply(dat.cells.sum.lst$K9m3$plate, function(p) plates.to.combine[[p]])


ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  facet_wrap(~platemerged) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

# can we fit? 

fit1 <- lm(formula = log2(totalcuts) ~ platemerged + cluster, data = dat.cells.sum.lst$K9m3)
fit2 <- lm(formula = log2(totalcuts) ~ plate + cluster, data = dat.cells.sum.lst$K9m3)

# no plate effects if we plot the prob
ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob / (1 - prob))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  facet_wrap(~plate) + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / K9") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob / (1 - prob))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / K9") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = log2(prob / (1 - prob)))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("log2(K36 / K9)") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Fit prob vs count -------------------------------------------------------


fit.prob <- lm(formula = prob / (1 - prob) ~ plate + cluster, data = dat.prob.avg)

dat.cells.sum.merged <- dat.cells.sum.lst %>% 
  bind_rows() %>%
  reshape2::dcast(., cluster ~ mark, value.var = "totalcuts", fun.aggregate = sum) %>%
  mutate(counts.ratio = K36 / (K9m3))
  # reshape2::dcast(., cell + cluster ~ mark, value.var = "totalcuts")

fit.counts <- lm(formula = log2(counts.ratio) ~ cluster, data = dat.cells.sum.merged)

ggplot(dat.cells.sum.merged, aes(x = cluster, y = log2(counts.ratio))) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dat.fits.prob <- data.frame(coefname = rownames(summary(fit.prob)$coefficients), 
                            summary(fit.prob)$coefficients, stringsAsFactors = FALSE) %>%
  filter(grepl(pattern = "^cluster", x = coefname))

# rename and reorder
dat.fits.prob <- dat.fits.prob %>%
  rowwise() %>%
  mutate(coefname = gsub(pattern = "^cluster", replacement = "", x = coefname))
coeffactors <- levels(dat.meta.filt$cluster)
coeffactors <- coeffactors[which(coeffactors != "Erythroid")]
dat.fits.prob$coefname <- factor(dat.fits.prob$coefname, levels = coeffactors)


ggplot(dat.fits.prob, aes(x = coefname, y = Estimate, ymin = Estimate - Std..Error, ymax = Estimate + Std..Error)) + 
  geom_point() + 
  geom_errorbar() + 
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  expand_limits(y = 0) + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


fit.k36.counts <- lm(formula = log2(totalcuts) ~ plate + cluster, data = dat.cells.sum.lst$K36)
fit.k9.counts <- lm(formula = log2(totalcuts) ~ plate + cluster, data = dat.cells.sum.lst$K9)

dat.fits.counts.k36 <- data.frame(coefname = rownames(summary(fit.k36.counts)$coefficients), 
                            summary(fit.k36.counts)$coefficients, stringsAsFactors = FALSE) %>%
  filter(grepl(pattern = "^cluster", x = coefname))

dat.fits.counts.k9 <- data.frame(coefname = rownames(summary(fit.k9.counts)$coefficients), 
                            summary(fit.k9.counts)$coefficients, stringsAsFactors = FALSE) %>%
  filter(grepl(pattern = "^cluster", x = coefname))
colnames(dat.fits.counts.k9) <- paste(colnames(dat.fits.counts.k9), "K9", sep = "_")


# Estimate the ratio of the two -------------------------------------------

dat.fits.counts.merged <- left_join(dat.fits.counts.k36, dat.fits.counts.k9, by = c("coefname" = "coefname_K9")) %>%
  rowwise() %>%
  mutate(Estimate.merged = Estimate - Estimate_K9,
         StdErr.merged = sqrt(Std..Error^2 + Std..Error_K9^2))

jfac <- 2
m.counts <- ggplot(dat.fits.counts.merged, aes(x = coefname, y = Estimate.merged)) + 
  geom_point() +
  geom_errorbar(mapping = aes(ymin = Estimate.merged - jfac * StdErr.merged, ymax = Estimate.merged + jfac * StdErr.merged)) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m.probs <- ggplot(dat.fits.prob, aes(x = coefname, y = Estimate, color = coefname)) + 
  geom_point() + 
  geom_errorbar(mapping = aes(ymin = Estimate - jfac * Std..Error, ymax = Estimate + jfac * Std..Error)) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_color_manual(values = cbPalette.coefficients)  + 
  expand_limits(y = 0) + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")


m.probs.raw <- ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = log2(prob / (1 - prob)))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("log2(K36 / K9)") + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



m.k36.counts <- ggplot(dat.cells.sum.lst$K36, aes(x = cluster, y = log10(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

m.k9.counts <- ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log10(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")



# Refit: but for K9me3 you don't know the cell types  ---------------------


dat.k9.fewercelltypes <- dat.cells.sum.lst$K9 %>%
  rowwise() %>%
  mutate(cluster = ifelse(as.character(cluster) %in% c("Erythroid", "WhiteBloodCells"), as.character(cluster), "NonBlood"))
dat.k9.fewercelltypes$cluster <- factor(dat.k9.fewercelltypes$cluster, levels = c("Erythroid", "WhiteBloodCells", "NonBlood"))

cbPalette.fewercelltypes <- cbPalette
cbPalette.fewercelltypes[[3]] <- "grey80"

fit.k9.counts.fewercelltypes <- lm(formula = log2(totalcuts) ~ plate + cluster, data = dat.k9.fewercelltypes)

estimate.nonblood <- summary(fit.k9.counts.fewercelltypes)$coefficients["clusterNonBlood", ]["Estimate"]
se.nonblood <- summary(fit.k9.counts.fewercelltypes)$coefficients["clusterNonBlood", ]["Std. Error"]

estimate.wbc <- summary(fit.k9.counts.fewercelltypes)$coefficients["clusterWhiteBloodCells", ]["Estimate"]
se.wbc <- summary(fit.k9.counts.fewercelltypes)$coefficients["clusterWhiteBloodCells", ]["Std. Error"]

dat.fits.counts.k9.fewer <- data.frame(Estimate_K9_nonblood = estimate.nonblood, 
                                       StdErr_K9_nonblood = se.nonblood, 
                                       Estimate_K9_wbc = estimate.wbc,
                                       StdErr_K9_wbc = se.wbc, 
                                       stringsAsFactors = FALSE)

dat.fits.counts.merged.fewercelltypes <- dat.fits.counts.k36 %>%
  rowwise() %>%
  mutate(Estimate.merged = ifelse(coefname == "clusterWhiteBloodCells", Estimate - estimate.wbc, Estimate - estimate.nonblood), 
         StdErr.merged = ifelse(coefname == "clusterWhiteBloodCells", sqrt(Std..Error^2 + se.wbc^2), sqrt(Std..Error^2 + se.nonblood^2)))


m.counts.fewercelltypes <- ggplot(dat.fits.counts.merged.fewercelltypes, aes(x = coefname, y = Estimate.merged)) + 
  geom_point() +
  geom_errorbar(mapping = aes(ymin = Estimate.merged - jfac * StdErr.merged, ymax = Estimate.merged + jfac * StdErr.merged)) +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


m.k9.fewercelltypes <- ggplot(dat.k9.fewercelltypes, aes(x = cluster, y = log10(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette.fewercelltypes) + 
  theme_bw() + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")


pdf(outpdf, useDingbats = FALSE)
  print(m.counts)
  print(m.probs)
  print(m.probs.raw)
  print(m.k36.counts)
  print(m.k9.counts)
  print(m.counts.fewercelltypes)
  print(m.k9.fewercelltypes)
  print(m.counts)
dev.off()


# Plot raw ----------------------------------------------------------------





