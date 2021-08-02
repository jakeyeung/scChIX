# Jake Yeung
# Date of Creation: 2021-07-21
# File: ~/projects/scChIX/analysis_scripts/3e-remove_blood_clean_cells_manual2.K27.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load mat  ---------------------------------------------------------------

jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

prefix <- "var_filtered_manual2"
prefix.out <- "var_filtered_manual2noblood"

hubprefix <- "/home/jyeung/hub_oudenaarden"

# inmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/var_filtered_manual2"
inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA")
dname <- paste0(prefix, "_", jstr)
indir <- file.path(inmain, dname)

outmain <- file.path("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2", prefix.out)
dir.create(outmain)
outdir <- file.path(outmain, jstr)
dir.create(outdir)

dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("celltyping_output_filt.", jmark, ".2021-07-20.rds")
  readRDS(file.path(indir, fname))
})

mats.lst <- lapply(jmarks, function(jmark){
  fname <- paste0("countmat_output_filt.", jmark, ".2021-07-20.rds")
  readRDS(file.path(indir, fname))
})

# bad clsts
jmark <- jmarks[[1]]
bad.clsts <- c("cluster10")

jmark <- jmarks[[2]]
bad.clsts <- c("cluster8")

jmark <- jmarks[[3]]
bad.clsts <- c("cluster8")

ggplot(dat.metas[[jmark]] %>% filter(cluster != bad.clsts), aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  ggtitle(jmark) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bad.clsts.lst <- list("cluster10", "cluster8", "cluster8")
names(bad.clsts.lst) <- names(jmarks)


# Write counts  ------------------------------------------------------------

bad.cells.lst <- lapply(jmarks, function(jmark){
  bad.cells <- subset(dat.metas[[jmark]], ! cluster %in% bad.clsts.lst[[jmark]])$cell
})

mats.filt.lst <- lapply(jmarks, function(jmark){
  jcells.keep <- bad.cells.lst[[jmark]]
  cols.keep <- colnames(mats.lst[[jmark]]) %in% jcells.keep
  mat.filt <- mats.lst[[jmark]][, cols.keep]
  print(paste("Removing blood cells:", ncol(mats.lst[[jmark]]) - ncol(mat.filt)))
  return(mat.filt)
})

metas.filt.lst <- lapply(jmarks, function(jmark){
  jcells.keep <- bad.cells.lst[[jmark]]
  meta.filt <- subset(dat.metas[[jmark]], cell %in% jcells.keep)
  return(meta.filt)
})

m.lst <- lapply(jmarks, function(jmark){
  ggplot(metas.filt.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

for (jmark in jmarks){
  print(jmark)
  outf.mat <- file.path(outdir, paste0("countmat_var_filt.", jmark, ".", Sys.Date(), ".rds"))
  outf.meta <- file.path(outdir, paste0("meta_var_filt.", jmark, ".", Sys.Date(), ".rds"))
  mat.filt.tmp <- mats.filt.lst[[jmark]]
  dat.meta.tmp <- metas.filt.lst[[jmark]]
  print(dim(mat.filt.tmp))
  print(dim(dat.meta.tmp))
  assertthat::assert_that(ncol(mat.filt.tmp) == nrow(dat.meta.tmp))
  saveRDS(mat.filt.tmp, file = outf.mat)
  saveRDS(dat.meta.tmp, file = outf.meta)
}
