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


inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata/metadata_cleaned.2021-10-04.rds"


# Clean up cluster names  -------------------------------------------------

dat.meta <- readRDS(inf.meta)

clstrs.order.orig <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "Stromal", "ConnectiveTissueProg")
names(clstrs.order.orig) <- clstrs.order.orig
clstrs.order <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "MesenchymalProgs", "Cardiomyocytes")

clsts.hash <- hash(clstrs.order.orig, clstrs.order)

dat.meta.cleaned <- dat.meta %>%
  rowwise() %>%
  mutate(cluster = scchicFuncs::AssignHash(x = as.character(cluster), jhash = clsts.hash, null.fill = as.character(cluster)))
dat.meta.cleaned$cluster <- factor(dat.meta.cleaned$cluster, levels = clstrs.order)


# Clean up cluster names in metadata for GEO  -----------------------------

inf.processed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/processed_data/scChIX_gastrulation_H3K36me3_and_H3K9me3_metadata.txt.gz"
dat.processed <- fread(inf.processed)

# update cluster
dat.processed$cluster <- sapply(dat.processed$cluster, function(x) clsts.hash[[x]]) 

# save new with new name 
table(dat.processed$cluster)



# Save  -------------------------------------------------------------------

outf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata/metadata_final_cluster_names.", Sys.Date(), ".txt")
fwrite(dat.meta.cleaned, file = outf.meta, sep = "\t")

outf.processed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/scChIX/processed_data/scChIX_gastrulation_H3K36me3_and_H3K9me3_metadata_final.txt"
# save processed output
fwrite(dat.processed, file = outf.processed, quote = TRUE, sep = "\t")

system(paste0("gzip ", outf.processed))
