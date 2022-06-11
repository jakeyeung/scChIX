# Jake Yeung
# Date of Creation: 2021-10-22
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/7-analyze_coords_dbl.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scChIX)
library(scchicFuncs)
library(ggforce)

hubprefix <- "/Users/yeung/hub_oudenaarden"

# Load dbl  ---------------------------------------------------------------

inf.output <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_K36-K9m3.RData")
load(inf.output, v=T)
fits.out <- act.repress.coord.lst


# Load meta ---------------------------------------------------------------

inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/metadata_flipped.rds"))
dat.meta <- readRDS(inf.meta)



# Wrangle -----------------------------------------------------------------


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


# Make 2d plot  -----------------------------------------------------------

m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)


# Annotate clusters -------------------------------------------------------

dat.dbl.annot <- subset(dat.meta, type == "dbl") %>%
  left_join(., coords.dbl) 

dat.dbl.summary <- dat.dbl.annot %>%
  group_by(louv.act, cluster) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louv.act) %>%
  mutate(nfrac = ncells / sum(ncells)) %>%
  arrange(desc(nfrac)) %>%
  group_by(louv.act) %>%
  filter(row_number() == 1)

print(dat.dbl.summary)

# rename 
jname <- c("ConnectiveTissueProg" = "MesenchymalProgs")

dat.dbl.summary <- dat.dbl.summary %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "ConnectiveTissueProg", "MesenchymalProgs", cluster)) %>%
  mutate(cluster = ifelse(cluster == "SchwannCellPrecusor", "SchwannCellPrecursor", cluster))

clst2celltype <- hash::hash(dat.dbl.summary$louv.act, dat.dbl.summary$cluster)


dat.dbl.summary.k9 <- dat.dbl.annot %>%
  group_by(louv.repress, cluster) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louv.repress) %>%
  mutate(nfrac = ncells / sum(ncells)) %>%
  arrange(desc(nfrac)) %>%
  group_by(louv.repress) %>%
  filter(row_number() == 1) %>%
  mutate(cluster = ifelse(cluster %in% c("Erythroid", "WhiteBloodCells"), cluster, "NonBlood"))


clst2k9 <- hash::hash(dat.dbl.summary.k9$louv.repress, dat.dbl.summary.k9$cluster)
  

# Replot  -----------------------------------------------------------------

coords.dbl.annot <- coords.dbl %>%
  rowwise() %>%
  mutate(cluster.k36 = clst2celltype[[louv.act]], 
         cluster.k9 = clst2k9[[louv.repress]]) %>%
  left_join(., dat.meta) %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "ConnectiveTissueProg", "MesenchymalProgs", cluster)) %>%
  mutate(cluster = ifelse(cluster == "SchwannCellPrecusor", "SchwannCellPrecursor", cluster))

# order them in a sane way
ctypes.k36 <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "Stromal", "MesenchymalProgs")
ctypes.k9 <- c("Erythroid", "WhiteBloodCells", "NonBlood")

# reorder
coords.dbl.annot$cluster.k36 <- factor(coords.dbl.annot$cluster.k36, levels = ctypes.k36)
coords.dbl.annot$cluster.k9 <- factor(coords.dbl.annot$cluster.k9, levels = ctypes.k9)
coords.dbl.annot$cluster <- factor(coords.dbl.annot$cluster, levels = ctypes.k36)

 
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597") 
clst2color <- hash::hash(ctypes.k36, cbPalette[1:length(ctypes.k36)])

coords.dbl.annot$colorcode <- sapply(coords.dbl.annot$cluster, function(x) clst2color[[as.character(x)]])

m.grid.reordered <- ggplot(coords.dbl.annot, aes(x = cluster.k36, y = cluster.k9, color = cluster)) +
  geom_point(alpha = 0.75, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.reordered)



# Write output ------------------------------------------------------------

outdir <- "/Users/yeung/data/dblchic/gastrulation"
outpdf <- file.path(outdir, paste0("H3K36me3_H3K9me3_downstream_plots/dbl_cell_assignment_to_cluster.", Sys.Date(), ".pdf"))
outrds <- file.path(outdir, paste0("H3K36me3_H3K9me3_downstream_objs/dbl_cell_assignment_to_cluster.", Sys.Date(), ".rds"))

pdf(outpdf, useDingbats = FALSE)
  print(m.grid.reordered)
dev.off()

# write metadata
saveRDS(coords.dbl.annot, file = outrds)



