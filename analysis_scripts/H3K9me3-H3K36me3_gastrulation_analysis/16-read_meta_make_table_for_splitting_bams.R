# Jake Yeung
# Date of Creation: 2022-07-14
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/16-read_meta_make_table_for_splitting_bams.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# split bams by Single or Double

jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks
jmarksall <- c("K36", "K9m3", "K36-K9m3"); names(jmarks) <- jmarks

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata"

dat.metas <- lapply(jmarks, function(jmark){
  fread(file.path(indir, paste0("metadata_cell_cluster_with_clustercol.", jmark, ".2021-12-03.txt")))
})

dat.metas.singles <- lapply(dat.metas, function(jdat){
  jdat.sub <- subset(jdat, type == "single", select = c(cell, type, cluster, stage))
  jdat.sub <- jdat.sub[!duplicated(jdat.sub), ]
  return(jdat.sub)
})

dat.metas.dbls <- lapply(dat.metas, function(jdat){
  subset(jdat, type == "dbl", select = c(cell, type, cluster, stage))
}) %>%
  bind_rows()
dat.metas.dbls <- dat.metas.dbls[!duplicated(dat.metas.dbls), ]

dat.metas.all <- c(dat.metas.singles, list("K36-K9m3" = dat.metas.dbls))

lapply(dat.metas.all, dim)

# First column cell, second column type  ----------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata"
for (jmark in jmarksall){
  fwrite(dat.metas.all[[jmark]], file.path(outdir, paste0("metadata_cell_experi.", jmark, ".", Sys.Date(), ".txt")), sep = "\t")
}


