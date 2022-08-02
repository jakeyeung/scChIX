# Jake Yeung
# Date of Creation: 2022-08-01
# File: ~/projects/scChIX/analysis_scripts/differentiation/1-setup_objs_for_notebook.R
# 

rm(list=ls())

set.seed(0)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)
library(igraph)
library(umap)

mark1="H3K4me1"
mark2="H3K36me3"
inflda="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_jupyterhub/lda_outputs_three_mat_GranuMacro_GrepGranulocyte_direct_filter/ldaAnalysis_dynamicgenes/lda_outputs.countmat_GranuMacro_DynamicGenes.H3K4me1-H3K36me3/ldaOut.countmat_GranuMacro_DynamicGenes.H3K4me1-H3K36me3.Robj"
lowess1="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_jupyterhub/outputs_for_glmpca/ptime_loessfits/Grep_Granulocytes_outlierfilt_H3K4me1_H3K36me3_H3K4me1-H3K36me3/ptime_lowess_fits_spline_linear.H3K4me1.2022-06-30.rds"
lowess2="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_jupyterhub/outputs_for_glmpca/ptime_loessfits/Grep_Granulocytes_outlierfilt_H3K4me1_H3K36me3_H3K4me1-H3K36me3/ptime_lowess_fits_spline_linear.H3K36me3.2022-06-30.rds"
w=0.77
inffits <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_jupyterhub/pseudotime_scchix_fits_fixedw_again_knn_all_args_Grep_Granulocytes_outlierfilt_multicore_gpu/H3K4me1_H3K36me3_H3K4me1-H3K36me3/scchix_fits_knn_25_H3K4me1-H3K36me3pseudotime_fixedw.0.77.2022-07-01.rds"
dat.fits.raw <- readRDS(inffits)


jmark1 <- "H3K4me1"
jmark2 <- "H3K36me3"
inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/from_jupyterhub/scchix_downstream_dbl_pseudotime_pretty_figures_KNN_GrepGranulocytes_outlierfilt/knn25/metadata_dat_merged_knn25.2022-07-02.txt"
dat.merge.rbind <- fread(inf.meta) %>%
  group_by(mark) %>%
  mutate(umap1.scale = 1 * scale(umap1, center = TRUE, scale = TRUE),
         umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
         umap1.shift = ifelse(mark == jmark1, umap1.scale - 3, umap1.scale + 3),
         umap2.scale = ifelse(mark == jmark1, -1 * umap2.scale, umap2.scale)) %>%
  arrange(daystr)

# clean up 
dat.merge.rbind.clean <- dat.merge.rbind %>%
  dplyr::select(c(cell, experi, mark, daystr, umap2.scale, umap1.shift)) %>%
  dplyr::rename(umap1 = umap1.shift, 
                umap2 = umap2.scale)

library(ggrastr)
m.joint <- ggplot(dat.merge.rbind.clean, aes(x = umap1, y = umap2, group = cell, color = daystr)) +
  geom_path(alpha = 0.05) +
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.5) +
  geom_point(alpha = 0.5) + 
  ggtitle("Single and double incubated cells together") +
  scale_color_viridis_d() +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.joint)


# Load LDA and do LDA  ----------------------------------------------------

load(inflda, v=T)


tm <- topicmodels::posterior(out.lda)
topics.mat <- tm$topics

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 50
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

umap.out <- umap(topics.mat, config = jsettings)



# Set up matdbl to get sum across neraest neighbors -----------------------

mat.dbl <- count.mat
nn <- 25
if (nn > 1){
  
  dbl.cells <- rownames(umap.out$knn$distances); names(dbl.cells) <- dbl.cells
  nearest.cells.lst <- lapply(dbl.cells, function(jcell){
    indx.vec <- umap.out$knn$indexes[jcell, 1:nn]
    nearest.cells <- rownames(umap.out$knn$distances)[indx.vec]
  })
  
  mat.dbl.knn <- lapply(dbl.cells, function(jcell){
    jcells.keep <- nearest.cells.lst[[jcell]]
    rowSums(mat.dbl[, jcells.keep])
  })
  mat.dbl <- do.call(cbind, mat.dbl.knn)
  
} else {
  print("nn = 1, no knn smoothing")
}

outdir <- "/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/data"

save(umap.out, file = "/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/data/MacDiffDblMatNearestNeighbors.RData")
# 
# # Save count mat ----------------------------------------------------------
# 
mat.dbl.knn <- mat.dbl
save(mat.dbl.knn, file = file.path(outdir, paste0("MacDiffDblMat_H3K4me1xH3K36me3.RData")))
# 
# # Save lowess -------------------------------------------------------------
# 
# lowess.fits.k4me1 <- readRDS(lowess1)
# lowess.fits.k36me3 <- readRDS(lowess2)
# 
# save(lowess.fits.k4me1, file = file.path(outdir, paste0("LowessFits_H3K4me1.RData")))
# save(lowess.fits.k36me3, file = file.path(outdir, paste0("LowessFits_H3K36me3.RData")))
# 
# 
# # Save fits raw -----------------------------------------------------------
# 
# save(dat.fits.raw, file = file.path(outdir, paste0("scChIXOutputs_H3K4me1xH3K36me3.RData")))
# 
# 
# # save meta ---------------------------------------------------------------
# 
# save(dat.merge.rbind.clean, file = file.path(outdir, paste0("MacDiffMetadata_H3K4me1xH3K36me3.RData")))
