# Jake Yeung
# Date of Creation: 2021-12-02
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/11-write_metadata_to_table_by_mark.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("K9m3", "K36"); names(jmarks) <- jmarks

# load metadata

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata/metadata_flipped.rds"

dat.meta <- readRDS(inf.meta)

dat.meta.lst<- split(dat.meta, dat.meta$mark)

# put just cell and cluster

dat.meta.filt.lst <- lapply(dat.meta.lst, function(jdat){
  subset(jdat, select = c(cell, cluster))
})


# Write output cell and cluster ------------------------------------------------------------

outdir <- "/Users/yeung/data/dblchic/gastrulation/metadata"
for(jmark in jmarks){
  outf.filt.tmp <- file.path(outdir, paste0("metadata_cell_cluster.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(dat.meta.filt.lst[[jmark]], file = outf.filt.tmp, sep = "\t")
}


# Write the rest all columns ----------------------------------------------------------

for(jmark in jmarks){
  outf.tmp <- file.path(outdir, paste0("metadata_cell_cluster_full.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(dat.meta.lst[[jmark]], file = outf.tmp, sep = "\t")
}




