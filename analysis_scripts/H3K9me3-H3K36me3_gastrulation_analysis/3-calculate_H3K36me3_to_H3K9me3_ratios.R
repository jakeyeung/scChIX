# Jake Yeung
# Date of Creation: 2021-09-30
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/3-calculate_H3K36me3_to_H3K9me3_ratios.R
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

pdf(file = file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/prob_summaries/dat_prob_avg_plots.", Sys.Date(), ".pdf")), useDingbats = FALSE)

# plot on the UMAP

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

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
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.prob.avg, aes(x = prob, fill = colorcode)) + 
  geom_density(alpha = 0.25) + 
  scale_fill_identity() + 
  # facet_wrap(~cluster, nrow = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# get avg for erythroid
eryth.avg <- median(subset(dat.prob.avg, cluster == "Erythroid")$prob)
ggplot(dat.prob.avg, aes(x = prob, fill = colorcode)) + 
  geom_density(alpha = 1) + 
  geom_vline(xintercept = eryth.avg, linetype = "dotted") + 
  scale_fill_identity() + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob)) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / (K36 + K9)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = log2(prob / (1 - prob)))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("log2(K36 / K9)") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.prob.avg, aes(x = cluster, fill = colorcode, y = prob / (1 - prob))) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab(label = "") + 
  ylab("K36 / K9") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()


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
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.cells.sum.lst$K36, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K36, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  facet_wrap(~plate) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = cluster, y = log2(totalcuts), fill = cluster)) + 
  geom_boxplot() + 
  facet_wrap(~plate) + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

ggplot(dat.cells.sum.lst$K9m3, aes(x = umap1, y = umap2, color = log2(totalcuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


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
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")

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
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

fit.prob <- lm(formula = prob / (1 - prob) ~ plate + cluster, data = dat.prob.avg)

# check by boxplot 

# plot fits with SE? 






# cnames.lst <- lapply(split(x = dat.meta.filt, f = dat.meta.filt$cluster), function(jdat){
#   jdat$cell
# })
# prob.mat.pbulk <- as.data.frame(scChIX::MeanAcrossClusters(count.mat = t(prob.mat), cnames.keep.lst = cnames.lst))


