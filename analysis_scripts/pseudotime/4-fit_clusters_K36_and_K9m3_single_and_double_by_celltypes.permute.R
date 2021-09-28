# Jake Yeung
# Date of Creation: 2021-08-24
# File: ~/projects/scChIX/analysis_scripts/pseudotime/4-fit_clusters_K36_and_K9m3_single_and_double_by_celltypes.R
# Load meta after unmixing and run DE to find neighborhood structures

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scChIX)


# Load raw counts (50kb genomewide) ---------------------------------------------------------

jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks

inf.lda.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj")
  assertthat::assert_that(file.exists(inf.lda.tmp))
  return(inf.lda.tmp)
})

out.objs <- lapply(inf.lda.lst, function(inf.lda){
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

count.mat.lst <- lapply(out.objs, function(jout){
  jout$count.mat
})


# Load meta ---------------------------------------------------------------

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_demux_cleaned_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/demux_cleaned_filtered_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3.2021-08-24.filt2.spread_7.single_and_dbl.txt"
assertthat::assert_that(file.exists(inf.meta))
dat.meta <- fread(inf.meta)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.meta, aes(x = umap1.shift, y = umap2.scale, color = cluster, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.01) +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.ctype <- dat.meta

cells.keep <- dat.ctype$cell

count.mat.filt.lst <- lapply(count.mat.lst, function(jcount){
  cols.keep <- colnames(jcount) %in% cells.keep
  jcount[, cols.keep]
})

# ggplot(dat.ctype, aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   facet_wrap(~type) +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# make Epithelial the reference celltype

dat.ctype <- dat.ctype %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "Epithelial", "aEpithelial", cluster))



# Fit each gene  ---------------------------------------------------------------


dat.annots.filt.mark <- dat.ctype
jname <- "manual2nocenterfilt2_K36_K9m3_K36-K9m3"
hubprefix <- "/home/jyeung/hub_oudenaarden"
ncores <- 16
jseed <- 123


outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/by_clusters", jname)
dir.create(outdir)
for (jmark in jmarks){

  outf <- file.path(outdir, paste0("glm_poisson_fits_output.clusters.", jname, ".", jmark, ".permute_seed_", jseed, ".RData"))

  count.mat <- count.mat.lst[[jmark]]

  set.seed(jseed)
  indx <- seq_len(ncol(count.mat))
  indx.permute <- sample(indx, size = length(indx), replace = FALSE)

  cnames.orig <- colnames(count.mat)
  cnames.permute <- cnames.orig[indx.permute]

  count.mat.permute <- count.mat
  colnames(count.mat.permute) <- cnames.permute

  cnames <- colnames(count.mat.permute)
  ncuts.cells.mark <- data.frame(cell = colnames(count.mat.permute), ncuts.total = colSums(count.mat.permute), stringsAsFactors = FALSE)

  jrow.names <- rownames(count.mat.permute)
  names(jrow.names) <- jrow.names

  print("fitting genes... permuted")
  system.time(
    jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
      jrow <- count.mat.permute[jrow.name, ]
      jout <- scChIX::FitGlmRowClusters.withse(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jrow.name, returnobj = FALSE, with.se = TRUE)
      return(jout)
    }, mc.cores = ncores)
  )
  save(jfits.lst, dat.annots.filt.mark, ncuts.cells.mark, count.mat.permute, dat.ctype, file = outf)
}







