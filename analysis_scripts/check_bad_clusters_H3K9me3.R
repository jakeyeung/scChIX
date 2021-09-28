# Jake Yeung
# Date of Creation: 2021-07-13
# File: ~/projects/scChIX/analysis_scripts/check_bad_clusters_H3K9me3.R
# Some plate-specific clusters?


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(scChIX)


# Load data  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks
outstr <- paste(jmarks, collapse = "_")

outmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt")
outdir <- file.path(outmain, outstr)
dir.create(outdir)

inmain <-  file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt")

dat.metaless.lst <- lapply(jmarks, function(jmark){
  inf.metaless <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt/celltyping_output_filt.", jmark, ".2021-07-13.rds"))
  readRDS(inf.metaless)
})

mat.metaless.lst <- lapply(jmarks, function(jmark){
  inf.mat.metaless <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt/countmat_output_filt.", jmark, ".2021-07-13.rds"))
  readRDS(inf.mat.metaless)
})

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.metaless.lst[[jmark]], aes(x = umap1, y = umap2, color = plate)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    # facet_wrap(~plate) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})


# Remove bad plates in H3K9me3 and rerun  ---------------------------------

bad.clusters.lst <- list("K36" = c(), "K9m3" = c("E10p5-B6C-K9m3-190409-2", "E9p5-CB6-K9m3-190409-1"), "K36-K9m3" = "E9p5-CB6-K36-K9m3-190409-2")

bad.plate.dbl <- "E9p5-CB6-K36-K9m3-190409-2"

dat.metaless.filt.k36 <- dat.metaless.lst$K36
dat.metaless.filt.k9m3 <- subset(dat.metaless.lst$K9m3, umap2 > -2)
dat.metaless.filt.dbl <- subset(dat.metaless.lst$`K36-K9m3`, plate != bad.clusters.dbl)

# dat.metaless.filt.lst <- lapply(dat.metaless.lst, function(jdat){
#   jdat.filt <- subset(jdat, !plate %in% unlist(bad.clusters.lst))
# })

dat.metaless.filt.lst <- list(dat.metaless.filt.k36, dat.metaless.filt.k9m3, dat.metaless.filt.dbl)
names(dat.metaless.filt.lst) <- jmarks

m.filt.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.metaless.filt.lst[[jmark]], aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(jmark) +
    facet_wrap(~plate) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

print(m.filt.lst)


# Write outputs -----------------------------------------------------------

mat.filt.lst <- lapply(jmarks, function(jmark){
  cells.keep <- dat.metaless.filt.lst[[jmark]]$cell
  jmat <- mat.metaless.lst[[jmark]]
  cols.keep <- colnames(jmat) %in% cells.keep
  jmat.filt <- jmat[, cols.keep]
  return(jmat.filt)
})


# write outputs and new metas?

jdate <- "2021-07-12"
for (jmark in jmarks){
  print(jmark)
  outmat <- file.path(outdir, paste0("count_tables.50000.", jmark, ".", jdate, ".rds"))
  outmeta <- file.path(outdir, paste0("meta_data.50000.", jmark, ".", jdate, ".rds"))
  mat.tmp <- mat.filt.lst[[jmark]]
  meta.tmp <- dat.metaless.filt.lst[[jmark]]
  print(dim(mat.tmp))
  print(dim(meta.tmp))
  saveRDS(mat.tmp, outmat)
  saveRDS(meta.tmp, outmeta)
}

