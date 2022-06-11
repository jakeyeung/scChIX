# Jake Yeung
# Date of Creation: 2022-01-02
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/15-plot_linked_UMAP_less_alpha.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load metadata -----------------------------------------------------------

outdir <- "/Users/yeung/data/dblchic/gastrulation/primetime_plots"

outpdf <- file.path(outdir, paste0("gastrulation_umap_final.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks

inf.meta.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata/metadata_cell_cluster_with_clustercol.", jmark, ".2021-12-03.txt")
  return(inf.tmp)
})

# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
dat.meta.long <- lapply(jmarks, function(jmark){
  dat <- fread(inf.meta.lst[[jmark]])
  if (jmark == "K9m3"){
    print("Shuffling K9me3 rows")
    rows <- sample(x = seq(nrow(dat)), size = nrow(dat), replace = FALSE)
    dat <- dat[rows, ]
  }
  return(dat)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(umap2.flip = ifelse(cluster == "WhiteBloodCells" & umap1.shift < 3.5, umap2.flip + 0.7, umap2.flip))

ggplot(dat.meta.long, aes(x = umap1.shift2, y = umap2.flip, color = clustercol, group = cell)) + 
  geom_path(alpha = 0.05) + 
  geom_point() + 
  scale_color_identity() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_minimal() + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# inf.rds <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/metadata_cleaned.2021-10-04.rds"
inf.rds <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata/metadata_flipped.rds"
dat.meta.rds <- readRDS(inf.rds) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(umap2.flip = ifelse(cluster == "WhiteBloodCells" & umap1.shift < 3.5, umap2.flip + 0.7, umap2.flip))
 
ggplot(dat.meta.rds, aes(x = umap1.shift2, y = umap2.flip, color = cluster, group = cell)) + 
  geom_path(alpha = 0.05) + 
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  # scale_color_manual(values = cbPalette) + 
  theme_minimal() + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

fwrite(dat.meta.long, paste0("metadata_cell_cluster_with_clustercol.bothmarks.primetime.", Sys.Date(), ".txt"), sep = "\t")