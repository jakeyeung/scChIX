# Jake Yeung
# Date of Creation: 2022-06-20
# File: ~/projects/scChIX/analysis_scripts/finals/get_processed_data.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load processed count tables ---------------------------------------------

jmarks <- paste(c("K36", "K9m3", "K36-K9m3")); names(jmarks) <- jmarks
jmarks.mac <- paste(c("H3K4me1", "H3K36me3", "H3K4me1-H3K36me3")); names(jmarks.mac) <- jmarks.mac

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/LDA_outputs_init"

infs <- lapply(jmarks, function(jmark){
  file.path(indir, paste0("ldaOut.countmat_var_filt.", jmark, ".Robj"))
})

outs <- lapply(infs, function(jinf){
  load(jinf, v=T)
  return(list(count.mat = count.mat, out.lda = out.lda))
})

# Save outputs ------------------------------------------------------------

count.mat.lst <- lapply(outs, function(jout){
  jout$count.mat
})

out.lda <- lapply(outs, function(jout){
  jout$out.lda
})


# Save outputs -----------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scchix_databank_GEO/processed_data"

writeouts <- lapply(jmarks, function(jmark){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_gastrulation.", jmark, ".Robj"))
  jmat <- count.mat.lst[[jmark]]
  saveRDS(object = jmat, file = outf)
})



# Get macrophage data  ----------------------------------------------------

indir.mac <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_cluster/snakemakes_round3/dynamicgenebody_GranuMacro_H3K4me1_H3K36me3_H3K4me1-H3K36me3/snakemake_outputs/LDA_outputs_init"

infs.mac <- lapply(jmarks.mac, function(jmark){
  file.path(indir.mac, paste0("ldaOut.countmat_var_filt.", jmark, ".Robj"))
})

outs.mac <- lapply(infs.mac, function(jinf){
  load(jinf, v=T)
  return(list(count.mat = count.mat, out.lda = out.lda))
})


writeouts.mac <- lapply(jmarks.mac, function(jmark){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_macrophagedifferentiation.", jmark, ".Robj"))
  jmat <- outs.mac[[jmark]]$count.mat
  saveRDS(object = jmat, file = outf)
})




