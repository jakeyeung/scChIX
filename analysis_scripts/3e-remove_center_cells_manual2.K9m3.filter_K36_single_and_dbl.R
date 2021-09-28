# Jake Yeung
# Date of Creation: 2021-08-02
# File: ~/projects/scChIX/analysis_scripts/3e-remove_center_cells_manual2.K9m3.filter_K36_single_and_dbl.R
# Remove mysterious neuron cluster


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load mat  ---------------------------------------------------------------

# jmarks <- c("K36", "K27", "K36-K27")
jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

prefix <- "var_filtered_manual2"
prefix.out <- "var_filtered_manual2nocenternoneu"

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


bad.clsts.lst <- list(c("cluster1"), c(), c())
names(bad.clsts.lst) <- names(jmarks)


# Write counts  ------------------------------------------------------------

good.cells.lst <- lapply(jmarks, function(jmark){
  good.cells <- subset(dat.metas[[jmark]], ! cluster %in% bad.clsts.lst[[jmark]])$cell
})


# Further define cells based after projecting cells  ----------------------

inf.meta.proj <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/celltyping_K36_first_try.2021-08-02.txt"
dat.meta.proj <- fread(inf.meta.proj)

bad.ctypes <- c("NeuralTubeNeuralProgs2", "NeuralTubeNeuralProgs3")
good.cells.noneu <- subset(dat.meta.proj, !celltype %in% bad.ctypes)$cell

good.cells.k36.before <- good.cells.lst$K36
print("K36 before")
print(length(good.cells.k36.before))
good.cells.k36.after <- good.cells.k36.before[good.cells.k36.before %in% good.cells.noneu]
print("K36 after")
print(length(good.cells.k36.after))
good.cells.lst$`K36` <- good.cells.k36.after

good.cells.dbl.before <- good.cells.lst$`K36-K9m3`
print("dbl before")
print(length(good.cells.dbl.before))
good.cells.dbl.after <- good.cells.dbl.before[good.cells.dbl.before %in% good.cells.noneu]
print("dbl after")
print(length(good.cells.dbl.after))
good.cells.lst$`K36-K9m3` <- good.cells.dbl.after


# Save outputs ------------------------------------------------------------




mats.filt.lst <- lapply(jmarks, function(jmark){
  jcells.keep <- good.cells.lst[[jmark]]
  cols.keep <- colnames(mats.lst[[jmark]]) %in% jcells.keep
  mat.filt <- mats.lst[[jmark]][, cols.keep]
  print(paste("Removing cells:", ncol(mats.lst[[jmark]]) - ncol(mat.filt)))
  return(mat.filt)
})

metas.filt.lst <- lapply(jmarks, function(jmark){
  jcells.keep <- good.cells.lst[[jmark]]
  meta.filt <- subset(dat.metas[[jmark]], cell %in% jcells.keep)
  return(meta.filt)
})

m.before.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.metas[[jmark]], aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
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
