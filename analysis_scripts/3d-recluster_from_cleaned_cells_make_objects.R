# Jake Yeung
# Date of Creation: 2021-07-10
# File: ~/projects/scChIX/analysis_scripts/3d-recluster_from_cleaned_cells_make_objects.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(scChIX)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# jsettings.louv <- umap.defaults
# jsettings.louv$n_neighbors <- 15
# jsettings.louv$min_dist <- 0.1
# jsettings.louv$random_state <- 123


# Load new LDA after cleaning  --------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

#
# jdate <- "2021-07-07"
# jsuffix <- "dbl_k36_cleaned"
# jsuffix <- "dbl_k36_k9m3_cleaned"
# jsuffix <- "dbl_k36_cleaned"-

jdate <- "2021-07-09"
jsuffix <- "dbl_cleaned"
jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks

for (jmark in jmarks){
  infrobj <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_", jsuffix, "/lda_outputs.countmat_output_filt.", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.countmat_output_filt.", jmark, ".", jdate, ".K-30.Robj"))
  load(infrobj, v=T)

  tm.result <- posterior(out.lda)

  dat.umap.long <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

  # annotate
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(plate, jsep = "-"),
           stage =  strsplit(cell, "-")[[1]][[1]],
           cluster = paste("cluster", louvain, sep = "")) %>%
    dplyr::select(-louvain)


  # Write tables  -----------------------------------------------------------

  # jdate <- Sys.Date()
  # write count tables, metadata, and LDA objects
  outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_cleaned_LDA2"
  dir.create(outmain)
  outdir <- file.path(outmain, jsuffix)
  dir.create(outdir)

  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste(jsuffix, jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  outmat <- file.path(outdir, paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outpdf <- file.path(outdir, paste0("plots.", jmark, ".", Sys.Date(), ".pdf"))

  saveRDS(count.mat, file = outmat)
  saveRDS(dat.umap.long, file = outmeta)
  saveRDS(out.lda, file = outlda)

  pdf(outpdf, useDingbats = FALSE)

  print(m)

  dev.off()



}






