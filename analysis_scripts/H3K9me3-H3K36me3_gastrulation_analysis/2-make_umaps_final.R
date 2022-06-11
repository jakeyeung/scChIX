# Jake Yeung
# Date of Creation: 2021-10-03
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/2-make_umaps_final.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load meta ---------------------------------------------------------------

hubprefix <- "/Users/yeung/hub_oudenaarden"

# inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/metadata_flipped.rds")
inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/dat_meta_ordered_colorcoded.rds")
outf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/metadata_cleaned/meatdata_cleaned.", Sys.Date(), ".rds"))

dat.meta <- readRDS(inf.meta) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(umap2.flip = ifelse(cluster == "WhiteBloodCells" & umap1.shift < 3.5, umap2.flip + 0.7, umap2.flip))

dat.meta.split <- split(dat.meta, dat.meta$mark) 

dat.meta.split$K9m3 <- dat.meta.split$K9m3 %>%
  arrange(desc(cluster))

dat.meta.merged <- dat.meta.split %>%
  bind_rows()

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette2 <- c("#ff9f7d", "#eb9d01", "#7fbedf")


clsts.remove <- c("Erythroid", "WhiteBloodCells")

inf.meta.color <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots/metadata_pseudotime_K36-K9me3.K9m3.2021-08-10.txt"))
assertthat::assert_that(file.exists(inf.meta.color))
dat.meta.color <- data.table::fread(inf.meta.color) %>%
  dplyr::select(cell, clustercol)

pdf(file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/umap_outputs/umap_outputs_H3K36me3_H3K9me3.", Sys.Date(), ".pdf")), useDingbats = FALSE)
ggplot(dat.meta, aes(x = umap1.shift2, y = umap2.flip, color = colorcode, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  theme_bw() + 
  scale_color_identity() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot just H3K9me3
dat.meta.k9 <- dat.meta %>% filter(mark == "K9m3" & !cluster %in% clsts.remove) 

ggplot(dat.meta.k9, aes(x = umap1.shift2, y = umap2.flip, color = stage, group = cell)) + 
  geom_point(alpha = 1) + 
  theme_minimal() + 
  scale_color_manual(values = c("#696969", "#9b870c", "#0072B2")) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1.shift2, y = umap2.flip, color = cluster, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


dat.meta.merged.cleaned <- subset(dat.meta.merged, select = c(cell, type, cluster, mark, stage, umap1.shift2, umap2.flip, colorcode)) %>%
  dplyr::rename(umap1 = umap1.shift2, 
                umap2 = umap2.flip)

saveRDS(dat.meta.merged.cleaned, file = outf.meta)  
  
  
