# Jake Yeung
# Date of Creation: 2022-07-19
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/18-read_textfiles_convert_to_RData.R
# Convert to RData and put inside scChIX package for vignettes 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("K36", "K9m3", "K36-K9m3"); names(jmarks) <- jmarks

# Load count mats ---------------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/scChIX/processed_data"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/data"
outf <- file.path(outdir, paste0("CountMatsGastrulationInputs.RData"))

countmats.gastru <- lapply(jmarks, function(jmark){
  ftmp <- file.path(indir, paste0("scChIX_gastrulation_", jmark, "_countmat.txt.gz"))
  indat <- fread(ftmp) %>%
    as.data.frame()
  rownames(indat) <- indat$V1
  indat$V1 <- NULL
  return(indat)
})

print(lapply(countmats.gastru, dim))



# Save obj ----------------------------------------------------------------

save(countmats.gastru, file = outf)
