# Jake Yeung
# Date of Creation: 2022-07-20
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/19-create_louvain_to_clustername_metadata.R
# For analyzing scChIX downstream, we need a louvain 2 cluster metadata


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata/metadata_final_cluster_names.2022-07-20.txt"
dat.meta <- fread(inf.meta)


# Load output -------------------------------------------------------------

inf.scchixoutput <- "/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/data/GastrulationScChIXOutputsK36K9m3.RData"

load(inf.scchixoutput, v=T)
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


# Make 2d plot  -----------------------------------------------------------

m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

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


dat.dbl.summary.k9 <- dat.dbl.annot %>%
  group_by(louv.repress, cluster) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louv.repress) %>%
  mutate(nfrac = ncells / sum(ncells)) %>%
  arrange(desc(nfrac)) %>%
  group_by(louv.repress) %>%
  filter(row_number() == 1) %>%
  mutate(cluster = ifelse(cluster %in% c("Erythroid", "WhiteBloodCells"), cluster, "NonBlood"))

print(dat.dbl.summary.k9)



# Finalize outputs --------------------------------------------------------

louvain.celltype.metadata <- list("K36" = dat.dbl.summary %>% select(louv.act, cluster), 
                                  "K9m3" = dat.dbl.summary.k9 %>% select(louv.repress, cluster))



# Save to scChIX package --------------------------------------------------

ctype.colcode.metadata <- subset(dat.meta, select = c(cluster, colorcode))
ctype.colcode.metadata <- ctype.colcode.metadata[!duplicated(ctype.colcode.metadata), ]

save(louvain.celltype.metadata, ctype.colcode.metadata, file = "/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/data/GastrulationLouvainCelltypeAnnotations.RData")
